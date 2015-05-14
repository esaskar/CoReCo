#!/usr/bin/env python
"""
Compute posterior probabilities for enzyme presence given a phylogenetic tree with BLAST and GTG scores for each species.
"""

import sys, math, datetime
import bisect

import common

from tree import newickToTree, findParentsAndChildren, getNewick, rootedTreeToStr
from bayesnet import BayesianNetwork
from digraph import Digraph

NODE_SUFFIX_BLAST = "_BLAST"
NODE_SUFFIX_GTG = "_GTG"

PARAM_P_ADDITION = "p-add"           # 0 -> 1
PARAM_P_DELETION = "p-del"           # 1 -> 0
DEFAULT_PARAMETERS = {PARAM_P_ADDITION : 0.1,
                      PARAM_P_DELETION : 0.1}

ROOT_NODE = "N1"

class MutationParameter:
    def __init__(self, p_add, p_del):
        self.p_add = p_add
        self.p_del = p_del

def read_internal_cpds(fn):
    f = open(fn)
    cpds = {}
    for u, v, p_add, p_del in common.read_lines(f):
        mp = MutationParameter(float(p_add), float(p_del))
        if u not in cpds:
            cpds[u] = {}
        cpds[u][v] = mp
    return cpds

def read_evidence(f):
    evidence = {}
    for node, reaction, blast, gtg in common.read_lines(f):
        if blast != "?":
            blast = float(blast)
        else:
            blast = 0.0
        if gtg != "?":
            gtg = float(gtg)
        else:
            gtg = 0.0
        if reaction not in evidence:
            evidence[reaction] = {}
        evidence[reaction][node] = {"BLAST" : blast, "GTG" : gtg}

    return evidence

def complete_evidence(evidence, observed_species):
    for r in evidence:
        for n in observed_species:
            if n not in evidence[r]:
                evidence[r][n] = {"BLAST" : 0, "GTG" : 0}

def read_parameters(ddir):
    global parameters
    f = open("%s/%s" % (ddir, common.PARAMETER_FILE))
    parameters = DEFAULT_PARAMETERS.copy()
    for s in f:
        if s.startswith("#"):
            continue
        vals = s.strip().split("\t")
        if vals[0] in DEFAULT_PARAMETERS:
            parameters[vals[0]] = vals[1]
        else:
            print "Ignoring parameter %s" % (vals[0])

def get_blast_node_name(u):
    return "%s%s" % (u, NODE_SUFFIX_BLAST)

def get_gtg_node_name(u):
    return "%s%s" % (u, NODE_SUFFIX_GTG)

def build_dag(T, use_tree = True, blast = True, gtg = True):
    parents, children = findParentsAndChildren(T, ROOT_NODE)
    dag = Digraph()
    observed_species = set()
    ancestors = set()
    if use_tree:
        for u in children:
            ancestors.add(u)
            for v in children[u]:
                dag.add_edge(u, v)
    for u in T:
        if u not in children:
            observed_species.add(u)
            if blast:
                dag.add_edge(u, get_blast_node_name(u))
            if gtg: 
                dag.add_edge(u, get_gtg_node_name(u))
    return dag, observed_species, ancestors

def set_internal_cpds(bnet, cpds):
    for u in bnet.dag.children:
        if len(bnet.dag.children[u]) > 0:
            if (len(bnet.dag.parents[u]) == 0):
                # root: uniform prior
                bnet.set_pd(u, [], {"0" : 0.5, "1" : 0.5})
            else:
                # species node
                parent = list(bnet.dag.parents[u])[0]

                mp = cpds[parent][u]
                
                bnet.set_pd(u, [parent], {"0" : {"0" : 1 - mp.p_add, 
                                                 "1" : mp.p_del}, 
                                          "1" : {"0" : mp.p_add, 
                                                 "1" : 1 - mp.p_del}})

def point_to_str(point):
    return "%f" % (point)

def parse_cpd(f):
    cpd = {}
    f.readline() # header
    for s in f:
        ix, x, density, p = s.strip().split()
        ix = int(ix.strip("\""))
        x = point_to_str(float(x))
        p = float(p)
        cpd[x] = p
    return cpd        

def merge_pos_and_neg_cpds(pos, neg):
    merged = {}
    xs = set()
    for x in pos:
        if x not in neg:
            raise Exception("Data point in positive but not in negative: %s" % (x))            
    for x in neg:
        if x not in pos:
            raise Exception("Data point in negative but not in positive %s" % (x))            

    for x in pos:
        merged[x] = {"0" : neg[x], "1" : pos[x]}
    return merged

def read_evidence_cpds(ddir):
    blastpos = parse_cpd(open("%s/blastpos" % (ddir)))
    blastneg = parse_cpd(open("%s/blastneg" % (ddir)))
    blast = merge_pos_and_neg_cpds(blastpos, blastneg)
    gtgpos = parse_cpd(open("%s/gtgpos" % (ddir)))
    gtgneg = parse_cpd(open("%s/gtgneg" % (ddir)))
    gtg = merge_pos_and_neg_cpds(gtgpos, gtgneg)
    return (blast, gtg)

def build_model(T, evidence_cpds, internal_cpds):
    dag, observed_species, ancestors = build_dag(T, use_tree = True, blast = True, gtg = True)
    bnet = BayesianNetwork()
    bnet.set_dag(dag)

    # set cpds for internal nodes
    set_internal_cpds(bnet, internal_cpds)

    # set cpds for evidence nodes
    blastcpd, gtgcpd = evidence_cpds
    for u in observed_species:
        bnet.set_pd(get_blast_node_name(u), [u], blastcpd)
        bnet.set_pd(get_gtg_node_name(u), [u], gtgcpd)

    return bnet, observed_species, ancestors

def build_blast_model(T, evidence_cpds, internal_cpds):
    dag, observed_species, ancestors = build_dag(T, use_tree = True, blast = True, gtg = False)
    bnet = BayesianNetwork()
    bnet.set_dag(dag)

    # set cpds for internal nodes
    set_internal_cpds(bnet, internal_cpds)

    # set cpds for evidence nodes
    blastcpd, gtgcpd = evidence_cpds
    for u in observed_species:
        bnet.set_pd(get_blast_node_name(u), [u], blastcpd)

    return bnet, observed_species, ancestors

def build_gtg_model(T, evidence_cpds, internal_cpds):
    dag, observed_species, ancestors = build_dag(T, use_tree = True, blast = False, gtg = True)
    bnet = BayesianNetwork()
    bnet.set_dag(dag)

    # set cpds for internal nodes
    set_internal_cpds(bnet, internal_cpds)

    # set cpds for evidence nodes
    blastcpd, gtgcpd = evidence_cpds
    for u in observed_species:
        bnet.set_pd(get_gtg_node_name(u), [u], gtgcpd)

    return bnet, observed_species, ancestors

def build_naive_model(T, cpds):
    dag, observed_species, ancestors = build_dag(T, use_tree = False, blast = True, gtg = True)
    bnet = BayesianNetwork()
    bnet.set_dag(dag)

    # set cpds for internal nodes
    for u in bnet.dag.children:
        if len(bnet.dag.children[u]) > 0 and len(bnet.dag.parents[u]) == 0:
            bnet.set_pd(u, [], {"0" : 0.5, "1" : 0.5})

    # set cpds for evidence nodes
    blastcpd, gtgcpd = cpds
    for u in observed_species:
        bnet.set_pd(get_blast_node_name(u), [u], blastcpd)
        bnet.set_pd(get_gtg_node_name(u), [u], gtgcpd)

    return bnet, observed_species, ancestors

def find_data_point(point, points):
    A = map(float, points)
    pos = bisect.bisect_left(A, float(point))
    if pos >= len(points):
        return points[-1]
    elif pos == 0:
        return points[0]
    else:
        if A[pos] - point < point - A[pos - 1]:
            return points[pos]
        else:
            return points[pos - 1]

def map_evidence(model, evidence, use_blast = True, use_gtg = True):
    mapped = {}
    for node in evidence:

        if node not in model.dag.parents:
            #print "Warning: evidence node %s not in tree" % (node)
            continue

        for source in evidence[node]:

            if source == "BLAST":
                if use_blast == False:
                    continue
                source_node = get_blast_node_name(node)
            elif source == "GTG":
                source_node = get_gtg_node_name(node)
                if use_gtg == False:
                    continue
            else:
                raise Exception("Invalid source in evidence: %s" % (source))

            points = model.pds[source_node].keys()
            points.sort(lambda x, y: cmp(float(x), float(y)))
            x = find_data_point(evidence[node][source], points)
            mapped[source_node] = {}
            for p in points:
                if p == x:
                    v = 1.0
                else:
                    v = 0.0
                mapped[source_node][p] = v
    return mapped

def compute_reaction_scores(model, naive_model, 
                            blast_model, gtg_model,
                            evidence, 
                            observed_species, ancestors, of, rrange = None):

    of.write("#Reaction Species Score NaiveScore ScoreLogRatio BlastTreeScore GTGTreeScore\n")
    c = 1
    keys = evidence.keys()
    keys.sort()
    if rrange != None:
        rmin, rmax = rrange
        keys = keys[rmin:rmax]
    for reaction in keys:
        sys.stderr.write("%d/%d: EC %s\n" % (c, len(evidence), reaction))
        mapped_evidence = map_evidence(model, evidence[reaction])
        ppd = model.compute_posterior(mapped_evidence)
        ppd_naive = naive_model.compute_posterior(mapped_evidence)

        mapped_evidence_blast = map_evidence(model, evidence[reaction], use_gtg = False)
        mapped_evidence_gtg = map_evidence(model, evidence[reaction], use_blast = False)        

        ppd_blast = blast_model.compute_posterior(mapped_evidence_blast)
        ppd_gtg = gtg_model.compute_posterior(mapped_evidence_gtg)

        k2 = ppd.keys()
        k2.sort()
        for node in k2:
            if node in observed_species or node in ancestors:
                score = common.Score()
                score.pscore = ppd[node]["1"]
                score.btscore = ppd_blast[node]["1"]
                score.gtscore = ppd_gtg[node]["1"]
                if score.btscore == None:
                    score.btscore = 0
                if score.gtscore == None:
                    score.gtscore = 0
                if node not in ppd_naive:
                    score.npscore = 0
                    score.logratio = 0
                else:
                    score.npscore = ppd_naive[node]["1"]
                    score.logratio = math.log(score.pscore / score.npscore) / math.log(2)
                of.write("%s\t%s\t%s\n" % (reaction, node, score))

        c += 1

def main(ddir, rrange = None):
    read_parameters(ddir)

    tf = open("%s/tree" % (ddir))
    t = newickToTree(tf.read())
    root = "N1"
    parent, children = findParentsAndChildren(t, root)
    treestr = rootedTreeToStr(parent, children, root)
    open("%s/tree.full" % (ddir), "w").write(treestr)

    evidence_cpds = read_evidence_cpds("%s/%s" % (ddir, common.CPD_DIR))
    internal_cpds = read_internal_cpds("%s/%s" % (ddir, common.TREE_CPD_FILE))
    
    ef = open("%s/evidence" % (ddir))
    evidence = read_evidence(ef)

    model, observed_species, ancestors = build_model(t, evidence_cpds, internal_cpds)
    naive_model, spec2, anc2 = build_naive_model(t, evidence_cpds)
    blast_model, spec3, anc3 = build_blast_model(t, evidence_cpds, internal_cpds)
    gtg_model, spec4, anc4 = build_gtg_model(t, evidence_cpds, internal_cpds)

    complete_evidence(evidence, observed_species)

    # compute posteriors
    of = open("%s/reaction-scores" % (ddir), "w")
    of.write("#Generated by \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))
    if rrange != None:
        of.write("#Reaction range: %s\n" % (rrange))
    compute_reaction_scores(model, naive_model, 
                            blast_model, gtg_model, evidence, 
                            observed_species, ancestors, of, rrange)

if __name__ == "__main__":
    if len(sys.argv) > 2:
        rrange = map(int, sys.argv[2].split(","))
    else:
        rrange = None
    main(sys.argv[1], rrange) # project dir


