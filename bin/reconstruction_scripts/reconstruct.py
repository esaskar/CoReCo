#!/usr/bin/env python

import sys, os, random, datetime, shutil

#import cvxopt
#import cvxopt.solvers
#from cvxopt import matrix, solvers
sys.path.append("../model_training_scripts/")

import atommap, common, network, kspyen
import visualize_atom_graph
#import network2sbml
import map_result_in_kegg

from priodict import priorityDictionary

from cost import rcost, min_cost, max_cost, cost_to_p, COST_EPSILON

import simple_graph

DEFAULT_REACTION_SCORE = 0.0

DEFAULT_ATOM_TYPES = ["C"]

# Add high-scoring reactions even if no gapfix found
ADD_FAILED_HIGHSCORE_REACTIONS = 1

# Add spontaenous reactions after reconstruction
ADD_SPONTANEOUS_REACTIONS = 1

MAX_RESULT_PATHWAYS = 25
MAX_VISITED_STATES = 500000
MAX_QUEUE_SIZE = 500
TERMINATE_ON_MAX_QUEUE = 0

DRAW_CANDIDATE_PATHWAYS = 0
DRAW_COMPLETE_PATHWAYS = 0
DRAW_PARTIAL_PATHWAYS = 0
DRAW_QUEUED_PATHWAYS = 0

# in seconds / target
MAX_SEARCH_TIME = 15 * 60

VERBOSE_NONE = 0
VERBOSE_SOME = 1
VERBOSE_ALL = 2

VERBOSE = 2

#NUM_PATHS = 50
#NUM_PATHS = 25
#NUM_PATHS = 25
NUM_PATHS = 5

#NUM_PATHS_IN_DEPTH = [200, 50]
#NUM_PATHS_IN_DEPTH = [1000, 50]
#K_PHASE_MULT = [1.0, 2.0]

REJECT_CYCLIC = 0

COMPUTE_YIELD = 1
REJECT_NONPOSITIVE_YIELD = 0
FBA_EPS = 1e-6

B = []

def readBoundaries(path): 
    f = open(path)   
    for s in f: 
	name, idreaction, model, lb, up = s.strip().split("\t")         
	bo = Bounds()
	bo.name = name
	bo.lb = lb
	bo.ub = up	   
	B.append(bo)
    f.close()

def searchBound(name):
    bound = B[1]	
    for b in B:
	if name in b.name:
	    bound = b	
            break
    return bound 
		
class Bounds:
    def __init__(self):
	self.name = 'bounds'
	self.lb = '-1000'
	self.up = '1000'

class Parameters:
    def __init__(self):
        pass

def message(level, msg):
    if VERBOSE >= level:
        sys.stdout.write("[%s] %s\n" % (datetime.datetime.now().strftime("%Y.%m.%d-%H:%M:%S"), msg))
        sys.stdout.flush()

def mols_for_atoms(A, amp):
    M = {}
    for a in A:
        mol, atom = amp.get_mol_and_atom(a)
        M[a] = mol
    return M

class SearchState:

    def __init__(self):
        self.R = set()
        self.targets = set()
        self.unsat = set()
        self.used_atoms = set()
        self.reached = set()
        self.cost = 0
        self.E = {}  # induced atom graph
        self.Erev = {} # reverse edges of self.E
        self.atom_paths = []
        self.tgt_yield = None

    def copy(self):
        s = SearchState()
        s.R = set(self.R)
        s.targets = set(self.targets)
        s.unsat = set(self.unsat)
        s.used_atoms = set(self.used_atoms)
        s.reached = set(self.reached)
        s.cost = self.cost
        s.atom_paths = list(self.atom_paths)
        s.tgt_yield = self.tgt_yield
        for u in self.E:
            s.E[u] = {}
            for v in self.E[u]:
                s.E[u][v] = self.E[u][v]
        for u in self.Erev:
            s.Erev[u] = {}
            for v in self.Erev[u]:
                s.Erev[u][v] = self.Erev[u][v]
        return s

    @classmethod
    def summary_header(cls):
        return "NumReactions TargetYield TotalCost AvgCost TotalLogP TotalP AvgLogP NumPaths Reactions"

    def average_cost(self):
        if len(self.R) > 0:
            return self.cost / len(self.R)
        else:
            return float("NAN")

    def total_log_prob(self):
        return self.cost - len(self.R) * COST_EPSILON

    def summarize(self):
        avgcost = self.average_cost()
        tlp = self.total_log_prob()
        tprob = cost_to_p(tlp)
        if len(self.R) > 0:
            avglp = tlp / len(self.R)
        else:
            avglp = float("NAN")
        rs = list(self.R)
        rs.sort()
        s = "%d\t%s\t%.3f\t%.3f\t%.3f\t%.5f\t%.3f\t%d\t%s" % (len(self.R), self.tgt_yield, self.cost, avgcost, tlp, tprob, avglp, len(self.atom_paths), ",".join(rs))
        return s

    def complete(self):
        return len(self.unsat) == 0

    def add_target(self, a):
        self.targets.add(a)
        self.add_unsat_atom(a)

    def add_unsat_atom(self, a):
        self.used_atoms.add(a)
        self.unsat.add(a)

    def add_used_atom(self, a):
        self.used_atoms.add(a)

    def add_reaction(self, r, amp, scores):
        self.add_path([r], amp, scores)

    def cost_estimate(self, atom_scores):
        if len(self.unsat) == 0:
            return 0
        estimate = float("inf")
        for u in self.unsat:
            if atom_scores[u] < estimate:
                estimate = atom_scores[u]
        return estimate

    def cost_estimate_indep(self, atom_scores):
        if len(self.unsat) == 0:
            return 0
        estimate = 0
	try:
            for u in self.unsat:
                estimate += atom_scores[u]
	except:
	    return 0 	
        return estimate

    def total_cost(self, atom_scores):
        return self.cost + self.cost_estimate(atom_scores)

    def update_cost(self, scores):
        self.cost = 0
        for r in self.R:
            self.cost += rcost(scores[r].pscore)

    def add_path(self, P, amp, scores):
        """Add reactions of P to the pathway and update pathway cost, unsatisfied and used atoms. Return True if updated pathway is cyclic and False otherwise."""
        self.R.update(P)
        self.update_cost(scores)
        cyclic = False
        used = set(self.used_atoms)
        edges = []

        # Determine the directed atom graph induced by the atom maps of reactions in P
        for r in P:
            re = amp.reactions[r]
            #print re, re.substrates, re.products
            for u in re.rev_maps:
                for v in re.rev_maps[u]:
                    if v not in self.E:
                        self.E[v] = {}
                    if u not in self.E:
                        self.E[u] = {}
                    if v not in self.Erev:
                        self.Erev[v] = {}
                    if u not in self.Erev:
                        self.Erev[u] = {}
                    self.E[v][u] = 1
                    self.Erev[u][v] = 1
                    edges.append((v, u))
                    used.add(u)
                    used.add(v)

        # is the induced directed atom graph cyclic?
        if REJECT_CYCLIC:
            #print "TARGETS", self.targets
            #print "Erev", self.Erev
            ends, Esub = simple_graph.backtrack(self.Erev, self.unsat)
            #print "Esub", Esub
            cyclic = simple_graph.is_cyclic(Esub)
            if len(Esub) > 0:
                tt = os.times()[0]
                simple_graph.quickdraw(self.E, "esub/Esub%s-E.svg" % (tt))
                simple_graph.quickdraw(Esub, "esub/Esub%s-%s.svg" % (tt, cyclic))
            #cyclic = simple_graph.is_cyclic(self.Erev)
        else:
            cyclic = None

        # determine reachable atoms in induced atom graph
        self.reached = simple_graph.find_reachable(self.E, self.reached)

        # determine unsatisfied atoms
        #self.unsat = simple_graph.backtrack(self.Erev, self.unsat, source_atoms).difference(source_atoms)
        #self.unsat = simple_graph.backtrack(self.Erev, self.unsat, source_atoms)
        #print "=== ADDING", P
        #print "=== UNSAT BEFORE:", self.unsat, mols_for_atoms(self.unsat, amp)
        #print "=== TARGETS:", self.targets, mols_for_atoms(self.targets, amp)
        #self.unsat = simple_graph.scc_backtrack(self.Erev, self.unsat)
        #self.unsat = simple_graph.scc_backtrack(self.Erev, self.targets)
        #print "=== UNSAT AFTER:", self.unsat, mols_for_atoms(self.unsat, amp)
        #self.unsat.difference_update(source_atoms)

        #for u in self.targets:
        #    print u, amp.get_mol_and_atom(u), u in self.Erev

        self.unsat = simple_graph.scc_backtrack_sources(self.Erev, self.targets, self.reached)
        #print "=== UNSAT - SOURCES:", self.unsat, mols_for_atoms(self.unsat, amp)
        self.used_atoms = used
        self.atom_paths.append(edges)

        #print "R", list(self.R)
        #ss = set(["R00243_2_rev", "R00258_1", "R00214_1_rev"])
        #if ss.issubset(self.R):
        #    sys.exit()

        return cyclic

    def __hash__(self):
        rr = list(self.R)
        rr.sort()
        uu = list(self.unsat)
        uu.sort()
        
        #print "---", hash((tuple(rr), tuple(uu))), rr, uu
        return hash((tuple(rr), tuple(uu)))

    def __eq__(self, other):
        if len(self.R) != len(other.R) or len(self.unsat) != len(other.unsat):
            return False
        for r in self.R:
            if r not in other.R:
                return False
        for u in self.unsat:
            if u not in other.unsat:
                return False
        return True

    def __str__(self):
        s = "Reactions:\n"
        k = list(self.R)
        k.sort()
        for r in k:
            s += "  %s\n" % (r)
        s += "Unsats:\n"
        k = list(self.unsat)
        k.sort()
        for u in k:
            s += "  %s\n" % (u)
        return s.rstrip("\n")

    def __repr__(self):
        return ""

def calculate_max_re_distances(max_mol_dists, amp):
    dists = {}
    for r in amp.reactions:
        re = amp.reactions[r]
        d = 0
        for mol in re.substrates:
            if mol in max_mol_dists and max_mol_dists[mol] > d:
                d = max_mol_dists[mol]
        bname = re.basename()
        if bname not in dists or dists[bname] > d:
            dists[bname] = d
    return dists

def estimate_reaction_costs(atom_scores, reactions, amp):
    """Reaction cost is the sum of substrate atom costs."""
    dists = {}
    for r in reactions:
        re = amp.reactions[r]
        d = 0.0
        for mol in re.substrates:
            for u in amp.atom_types[mol]:
                if u in atom_scores:
                    d += atom_scores[u]
        bname = re.basename()
        if bname not in dists or dists[bname] > d:
            dists[bname] = d
    return dists

def write_parameters(params):
    o = open("%s/%s" % (params.odir, common.PARAMETER_FILE), "w")
    o.write("Threshold\t%s\n" % (params.threshold))
    o.write("PathRejectThreshold\t%s\n" % (params.path_reject_threshold))
    o.write("CostEpsilon\t%s\n" % (COST_EPSILON))
    o.write("MinCost\t%s\n" % (min_cost()))
    o.write("MaxCost\t%s\n" % (max_cost()))
    o.write("DefaultReactionScore\t%s\n" % (DEFAULT_REACTION_SCORE))
    o.write("MaxResultPathways\t%s\n" % (MAX_RESULT_PATHWAYS))
    o.write("MaxVisitedStates\t%s\n" % (MAX_VISITED_STATES))
    o.write("MaxQueueSize\t%s\n" % (MAX_QUEUE_SIZE))
    o.write("MaxSearchTime\t%s\n" % (MAX_SEARCH_TIME))
    o.write("NumPaths\t%s\n" % (NUM_PATHS))
    o.write("RejectCyclic\t%s\n" % (REJECT_CYCLIC))
    o.write("ComputeYield\t%s\n" % (COMPUTE_YIELD))
    o.write("RejectNonpositiveYield\t%s\n" % (REJECT_NONPOSITIVE_YIELD))
    o.write("FBAEpsilon\t%s\n" % (FBA_EPS))
    o.close()

def find_aux_sources(R, amp, sources):
    res = set()
    for r in R:
        re = amp.reactions[r]
        res.add(re)
    agi = network.AtomGraph(res, lambda x: 1) # atom graph induced by R
    ratoms = simple_graph.find_reachable(agi.E, sources)

def read_reaction_scores(fn, amp, r2ec, threshold):
    f = open(fn)
    scores = {}
    ecscores = {}
    for s in f:
        vals = s.strip().split("\t")
        rid, ec = vals[0:2]
        ec = ec.split(",")
        rscore = common.Score("\t".join(vals[2:]))
        rscore.ec = ec
        if rid in amp.base_reactions:
            for rid2 in amp.base_reactions[rid]:
                scores[rid2] = rscore
        else:
            #print "Scored reaction %s not in atom maps" % (rid)
            ecscores[",".join(ec)] = rscore

    for r in amp.reactions:
        if r not in scores:
            rscore = common.Score()
            if r not in r2ec or r2ec[r] == ["?"]:
		# Spontaneous reactions
                # No EC association known for this reaction
                rscore.pscore = threshold - 1e-6  # this should fix the problem of reactions appearing in models without gene evidence with score==threshold
            else:
                # EC association known but EC not scored (EC does not appear in reaction-scores.full)
                # does not appear in ec_files.txt (not in swiss-prot training data)
                rscore.pscore = DEFAULT_REACTION_SCORE
            rscore.ec = ["?"] # TODO
            rscore.npscore = rscore.btscore = rscore.gtscore = rscore.bscore = rscore.gscore = 0
            scores[r] = rscore
    return scores, ecscores

def propagate_reaction_scores(scores, equivalences):
    """Sets reaction score S(r) = max(S(r')) where reaction r and r' have EC numbers that belong to the same equivalence class.

    This is to alleviate problems with mismatching EC annotations between KEGG and UniProt. """

    #for r in scores:
    #    rs = scores[r]
    pass

def compute_atom_score_heuristic(sources = None, source_atoms = None, params = None):
    assert(sources == None or source_atoms == None and sources != source_atoms)
    assert(params != None)

    if sources != None:
        source_atoms = get_mol_atoms(sources, params.amp)

    atom_scores = {}
    preceding_reactions = {}
    Q = priorityDictionary()
    k = params.amp.atom_types.keys()
    k.sort()
    for index in source_atoms:
        atom_scores[index] = 0
        preceding_reactions[index] = set()
        Q[index] = 0
    visited = set(source_atoms)
    for q in Q:
        mol, atom = params.amp.get_mol_and_atom(q)
        atom_scores[q] = Q[q]
        visited.add(q)
        if mol in params.net.G:
            for r in params.net.G[mol]:
                if r in preceding_reactions[q]:
                    cost = 0
                else:
                    cost = rcost(params.scores[r].pscore)
                reaction = params.amp.reactions[r]
                #print r, q, mol
                if q in reaction.maps:
                    for mapped_atom in reaction.maps[q]:
                        new_cost = Q[q] + cost
                        if mapped_atom not in visited and (mapped_atom not in Q or new_cost < Q[mapped_atom]):
                            Q[mapped_atom] = new_cost
                            preceding_reactions[mapped_atom] = set([r]).union(preceding_reactions[q])

    return visited, atom_scores, preceding_reactions

def write_atom_score_distribution(atom_scores):
    o = open("atom-scores.txt", "w")
    for u in atom_scores:
        o.write("%s\n" % (atom_scores[u]))
    o.close()

def generate_scores(reactions):
    base_scores = {}
    scores = {}
    import random
    
    rnd = random.Random(0)
    for r in reactions:
        base = r.split("_")[0]
        if base not in base_scores:
            #base_scores[base] = rnd.uniform(0, 1)
            base_scores[base] = 0.5

        scores[r] = base_scores[base]
    return scores

def precalc_colors():
    global REACTION_COLORS
    REACTION_COLORS = []
    n = 50
    for i in range(n):
        r = 1.0 * i / n
        g = 1
        b = 0
        REACTION_COLORS.append((r, g, b))
    n = 200
    for i in range(n):
        r = 1
        g = 1.0 * (n - i) / n
        b = 0
        REACTION_COLORS.append((r, g, b))

def cost_to_color(x):
    global REACTION_COLORS
    minv = min_cost()
    maxv = max_cost()
    sx = (x - minv) / (maxv - minv)
    ix = int(sx * len(REACTION_COLORS))
    if ix < 0:
        ix = 0
    if ix >= len(REACTION_COLORS):
        ix = len(REACTION_COLORS) - 1
   
    r, g, b = REACTION_COLORS[ix]

    s = "#%.2x%.2x%.2x" % (r * 255, g * 255, b * 255)
    return s

PATH_COLORS = ["red", "green", "blue", "orange", "seagreen", "purple", "cadetblue", "saddlebrown", "slategrey"]

def draw_state(agv, target_reactions, target_molecules, amp, rstate, sources, ofn, scores, reachable_atoms, draw_edge_labels = False, mol_names = {}, atom_scores = {}):

    rcolors = {}
    rfillcolors = {}
    for r in target_reactions:
        rcolors[r] = "orange"
    acolors = {}
    mcolors = {}
    mfillcolors = {}
    all_mols = set()

    atom_edge_colors = {}
    c = 0
    for p in rstate.atom_paths:
        for u, v in p:
            if u not in atom_edge_colors:
                atom_edge_colors[u] = {}
            atom_edge_colors[u][v] = PATH_COLORS[c]
        c = (c + 1) % len(PATH_COLORS)

    for u in rstate.used_atoms:
        acolors[u] = "blue"

    for u in rstate.reached:
        acolors[u] = "purple"

    for mol in sources:
        all_mols.add(mol)
        mfillcolors[mol] = "palegreen"
        if mol in amp.atom_types:
            for u in amp.atom_types[mol]:
                acolors[amp.get_atom_index(mol, u)] = "palegreen"

    for mol in target_molecules:
        all_mols.add(mol)
        mfillcolors[mol] = "orange"
        if mol in amp.atom_types:
            for u in amp.atom_types[mol]:
                acolors[amp.get_atom_index(mol, u)] = "orange"

    for u in rstate.unsat:
        acolors[u] = "red"
        mol, atom = amp.get_mol_and_atom(u)
        all_mols.add(mol)
        mfillcolors[mol] = "red"

    atoms_to_draw = set()
    path_cost = rstate.cost
    cost_est = rstate.cost_estimate(atom_scores)
    rlabels = {}
    for rr in rstate.R:

        scorelabel = "p:%.2f,n:%.2f\\nbt:%.2f,gt:%.2f\\nB:%.2f,G:%.2f" % (scores[rr].pscore,
                                           scores[rr].npscore,
                                           scores[rr].btscore,
                                           scores[rr].gtscore,
                                           scores[rr].bscore,
                                           scores[rr].gscore)

        re = amp.reactions[rr]
        cost = rcost(scores[rr].pscore)
        rfillcolors[rr] = cost_to_color(cost)
        #total_cost += cost
        eclabel = ",".join(scores[rr].ec[:min(len(scores[rr].ec), 3)])
        rlabels[rr] = "%s (%s)\\n%s\\n%s"  % (rr, eclabel, cost, scorelabel)
        mols = set(re.substrates)
        mols.update(re.products)
        all_mols.update(mols)
        for mol in mols:
            for a in amp.atom_types[mol]:
                u = amp.get_atom_index(mol, a)
                if u in reachable_atoms:
                    atoms_to_draw.add(u)

    mol_scores = {}
    for u in atoms_to_draw:
        mol, atom = amp.get_mol_and_atom(u)
        if u in atom_scores:
            sc = atom_scores[u]
        else:
            sc = float("inf")
        if mol not in mol_scores or mol_scores[mol] > sc:
            mol_scores[mol] = sc

    mlabels = {}
    for mol in all_mols:
        if mol in mol_names:
            mlabels[mol] = "%s\\n%s" % (mol, mol_names[mol])
        else:
            mlabels[mol] = mol
        if mol in mol_scores:
            mlabels[mol] += "\\n%.2f" % (mol_scores[mol])

    glabel = "%d reactions, %d paths, %d unsats, cost=%.3f, est=%.3f, total=%.3f, yield=%s" % (len(rstate.R), len(rstate.atom_paths), len(rstate.unsat), path_cost, cost_est, path_cost + cost_est, rstate.tgt_yield)

    agv.draw_metabolic_network(rstate.R, "%s_net" % (ofn), rcolors, rfillcolors, rlabels, mcolors, mfillcolors, mlabels, "svg", glabel)
    agv.draw_atom_graph(rstate.R, atoms_to_draw, "%s_ag" % (ofn), rcolors, acolors, "svg", draw_edge_labels, glabel, atom_edge_colors)
    #agv.draw_metabolic_network(rstate.R, "%s_net" % (ofn), rcolors, rfillcolors, rlabels, mcolors, mfillcolors, mlabels, "png", glabel)
    #agv.draw_atom_graph(rstate.R, atoms_to_draw, "%s_ag" % (ofn), rcolors, acolors, "png", draw_edge_labels, glabel, atom_edge_colors)

def visualize_path(target_reactions, target_molecules, i, path, sources, scores, amp, reachable_atoms, mol_names, atom_scores, odir, cdir):
    agv = visualize_atom_graph.AtomGraphVisualizer(amp, cdir)
    if len(target_reactions) > 0:
        ddir = "%s/%s/%s" % (odir, common.PATHWAYS_DIR, "-".join(target_reactions))
    elif len(target_molecules) > 0:
        ddir = "%s/%s/%s" % (odir, common.PATHWAYS_DIR, "-".join(target_molecules))
    else:
        assert(0) # both target_r and target_m empty
        
    try: os.mkdir("%s/%s" % (odir, common.PATHWAYS_DIR)) 
    except: pass
    try: os.mkdir(ddir) 
    except: pass

    draw_state(agv, target_reactions, target_molecules, amp, path, sources, "%s/%s" % (ddir, i), scores, reachable_atoms, draw_edge_labels = 0, mol_names = mol_names, atom_scores = atom_scores)

def get_mol_atoms(mols, amp):
    atoms = set()
    for mol in mols:
        if mol in amp.atom_types:
            for u in amp.atom_types[mol]:
                a = amp.get_atom_index(mol, u)
                atoms.add(a)
    return atoms

def reaction_id_map_dir(rname):
    """Reaction name -> (id, map-id, is-forward-reaction)"""
    vals = rname.split("_")
    if len(vals) == 2:
        return vals[0], int(vals[1]), True
    else:
        return vals[0], int(vals[1]), False

SRC_NODE = "src"
TGT_NODE = "tgt"

def contains_both_directions(p):
    bases = {}
    for r in p:
        rid, mid, fwd = reaction_id_map_dir(r)
        if rid not in bases:
            bases[rid] = fwd
        elif fwd != bases[rid]:
            return True
    return False

def find_atom_paths(k, src_atoms, tgt_atoms, params):
    """Find k shortest atom paths.

k -- #shortest paths to find
src_atoms -- source nodes
tgt_atoms -- target nodes
    """

    orig_edge_costs = {}

    # atom graph set up:
    # add target edges and collect edges outgoing from target nodes
    for u in tgt_atoms:
        #assert(u in params.reachable_atoms)
	try:
            if u not in orig_edge_costs:
                 orig_edge_costs[u] = {}
            for v in params.ag.E[u]:
                 orig_edge_costs[u][v] = params.ag.E[u][v]
            params.ag.add_edge(u, TGT_NODE, 0, [])
	except: pass	   

    # collect incoming edges to source and forbidden nodes
    # NOTE: forbidden_atoms is currently always empty
    src_forb = src_atoms.union(params.forbidden_atoms)
    for u in src_forb:
        if u not in tgt_atoms:
            if u in params.ag.Erev:
                # for each incoming edge to u
                for v in params.ag.Erev[u]:
                    if v not in orig_edge_costs:
                        orig_edge_costs[v] = {}
                    orig_edge_costs[v][u] = params.ag.E[v][u]

    # add source edges
    for u in src_atoms:
        params.ag.add_edge(SRC_NODE, u, 0, [])

    # delete collected, forbidden edges
    for u in orig_edge_costs:
        for v in orig_edge_costs[u]:
            del params.ag.E[u][v]

    paths = kspyen.kspSimpleYen(params.ag, k, SRC_NODE, TGT_NODE)

    # revert modifications to the atom graph
    for u in tgt_atoms:
        params.ag.delete_edge(u, TGT_NODE)
    for u in src_atoms:
        params.ag.delete_edge(SRC_NODE, u)

    for u in orig_edge_costs:
        for v in orig_edge_costs[u]:
            params.ag.E[u][v] = orig_edge_costs[u][v]

    if [] in paths:
        paths.remove([])

    noncyclics = []
    #tt = os.times()[0]
    for pi, p in enumerate(paths):
        atoms = set()
        #E = {}
        cyclic = False
        for i, u in enumerate(p[0:-1]):
            if u in atoms:
                cyclic = True
                break
            atoms.add(u)
            #v = p[i + 1]
            #if u not in E:
            #    E[u] = {}
            #E[u][v] = 1

        if not cyclic:
            noncyclics.append(p)

        #simple_graph.quickdraw(E, "esub/%s-path%d-%s.svg" % (tt, pi, cyclic))

    return noncyclics #paths


def collect_reaction_sets(paths, tgt_atoms, params):
    reaction_sets = set()
    for pi, p in enumerate(paths):
        # p is an atom path SRC_NODE -> p_1 -> p_2 -> ... -> p_k -> TGTNODE
        adds = []
        for i in range(1, len(p) - 2):
            u = p[i]
            if u in tgt_atoms:
                # reject the pathway if targets found also in the middle
                adds = []
                break
            v = p[i + 1]
            bestr = None
            bestc = float("inf")
            for cr in params.ag.edge_reactions[u][v]:
                ccost = rcost(params.scores[cr].pscore)
                if ccost < bestc:
                    bestr = cr
                    bestc = ccost
            adds.insert(0, bestr) 	    
	 
        if contains_both_directions(adds):
            #print "        REJECT - DIRECTIONALITY"#, adds
            continue
	
        reaction_sets.add(tuple(adds))

    return reaction_sets

# from cvxopt import matrix, solvers  
# c = matrix([-4., -5.])  
# G = matrix([[2., 1., -1., 0.], [1., 2., 0., -1.]])  
# h = matrix([3., 3., 0., 0.])  
# sol = solvers.lp(c, G, h)  
# print "\nx = \n\n", sol['x']  

# c = matrix([0., 0., 0., 0., 0., 0., -1.])
# G = matrix([[-1., 0., 0., 0., 0., 0., 0., 1., 0.],
#             [0., -1., 0., 0., 0., 0., 0., 0., 1.], 
#             [0., 0., -1., 0., 0., 0., 0., 0., 0.], 
#             [0., 0., 0., -1., 0., 0., 0., 0., 0.], 
#             [0., 0., 0., 0., -1., 0., 0., 0., 0.], 
#             [0., 0., 0., 0., 0., -1., 0., 0., 0.], 
#             [0., 0., 0., 0., 0., 0., -1., 0., 0.]])
# h = matrix([[0., 0., 0., 0., 0., 0., 0., 1., 1.]])
# A = matrix([[1., 0., 0., 0., 0., 0.],
#             [0., 1., 0., 0., 0., 0.],
#             [-1., 0., 1., 0., 0., 0.],
#             [0., -1., 0., 1., 0., 0.],
#             [0., 0., -1., -1., 1., 0.],
#             [0., 0., 0., 0., -1., 1.],
#             [0., 0., 0., 0., 0., -1.]])
# b = matrix([0. for x in range(6)])

# print c
# print G
# print h
# print A
# print b

# sol = solvers.lp(c, G, h, A, b)  

next_mol_index = 0
def get_index(idd, D):
    global next_mol_index
    if idd not in D:
        D[idd] = next_mol_index
        next_mol_index += 1
    return D[idd]

def rev_index(D):
    Dr = {}
    for k in D:
        Dr[D[k]] = k
    return Dr

def print_matrix(X):
    o = sys.stdout
    for row in range(X.size[0]):
        for col in range(X.size[1]):
            o.write("% .2f " % (X[row,col]))
        o.write("\n")

PRINT_FBA = 0
            
def compute_yield(state, target_molecules, sources, amp):
    # balanced metabolite <-> contain at least 1 reached atom that is not source
    balanced_mols = set()
    for u in state.reached:
        if u not in sources:
            mol, atom = amp.get_mol_and_atom(u)
            balanced_mols.add(mol)

    # source metabolite <-> contain at least 1 source atom
    source_mols = set()
    for u in sources:
        mol, atom = amp.get_mol_and_atom(u)
        source_mols.add(mol)

    all_mols = set(balanced_mols).union(source_mols)

    N = {} # stoichiometric matrix
    objv = {} # objective vector

    r2ix = {}
    m2ix = {}
    global next_mol_index
    next_mol_index = 0
    next_r_index = 0
    for r in state.R:
        N[next_r_index] = {}
        r2ix[r] = next_r_index
        #if PRINT_FBA:
        #    print r, "->", next_r_index
        re = amp.reactions[r]
        all_mols.update(re.substrates)
        all_mols.update(re.products)
        #print re.sub_coeff
        #print re.pro_coeff
        for mol in re.substrates:
            mol_ix = get_index(mol, m2ix)
            #N[next_r_index][mol_ix] = -1.0
            N[next_r_index][mol_ix] = -re.sub_coeff[mol]
        for mol in re.products:
            mol_ix = get_index(mol, m2ix)
            #N[next_r_index][mol_ix] = 1.0
            N[next_r_index][mol_ix] = re.pro_coeff[mol]
        next_r_index += 1

    # add metabolite sink reactions and set up objective vector
    for mol in all_mols:
        mol_ix = get_index(mol, m2ix)
        #if PRINT_FBA:
        #    print "SINK", mol, mol_ix, "->", next_r_index
        N[next_r_index] = {}
        N[next_r_index][mol_ix] = -1.0
        if mol in target_molecules:
            objv[next_r_index] = -1.0
        next_r_index += 1

    not_balanced = all_mols.difference(balanced_mols)
    # if PRINT_FBA:
    #     print "NOT_BALANCED", not_balanced
    # add metabolite source reactions
    intakes = {}
    for mol in not_balanced:
        mol_ix = get_index(mol, m2ix)
        #if PRINT_FBA:
        #    print "SOURCE", mol, mol_ix, "->", next_r_index
        N[next_r_index] = {}
        N[next_r_index][mol_ix] = 1.0
        intakes[next_r_index] = mol_ix
        next_r_index += 1

    ix2r = rev_index(r2ix)
    ix2m = rev_index(m2ix)
    num_r = next_r_index
    num_m = next_mol_index

    #print ix2m

    #print "==== N ====\n", N
    #print "MAXR", num_r
    #print "MAXM", num_m

    # set up linear program for FBA

    # objective            -> c
    # reaction rates >= 0  -> G + h
    # steady-state         -> A + b
    # intake 0 <= x <= 1

    # c = [1 x r]
    c = cvxopt.matrix([0.0 for x in range(num_r)])
    for i in objv:
        c[i] = objv[i]

    # G = [r x r]
    #G = cvxopt.spmatrix(-1.0, range(num_r + len(intakes)), range(num_r))
    G = cvxopt.spmatrix(0.0, [], [], (num_r + len(intakes), num_r))
    for i in range(num_r):
        G[i, i] = -1.0
    # h = [1 x r]
    h = cvxopt.matrix([0.0 for x in range(num_r + len(intakes))])

    inix = num_r
    for ri in intakes:
        mi = intakes[ri]
        G[inix, ri] = 1
        h[inix] = 1
        inix += 1

    # A = [m x r]
    #num_used_mols
    A = cvxopt.spmatrix(0.0, [num_m - 1], [num_r - 1])
    for ri in N:
        for mi in N[ri]:
            A[mi, ri] = N[ri][mi]

    # b = [1 x m]
    b = cvxopt.matrix([0.0 for x in range(num_m)])

    #print "c",c
    #print "G"
    #print_matrix(G)
    #print "h",h
    #print "A"
    #print_matrix(A)
    #print "b",b

    #print "c", c.size
    #print "G", G.size
    #print "h", h.size
    #print "A", A.size
    #print "b", b.size

    cvxopt.solvers.options["show_progress"] = False
    sol = cvxopt.solvers.lp(c, G, h, A, b)
    if sol["status"] != "optimal":
        print "ERROR - result not optimal"
        return 0
    
    res = -sol["primal objective"]
    if res < FBA_EPS:
        res = 0
    return res


def filter_paths(paths, amp):
    filtered = []
    for pi, p in enumerate(paths):
        visited_mols = set()
        cyclic = False
        # p is an atom path SRC_NODE -> p_1 -> p_2 -> ... -> p_k -> TGTNODE
        for i in range(1, len(p) - 2):
            u = p[i]
            mol, atom = amp.get_mol_and_atom(u)
            if mol in visited_mols:
                cyclic = True
                break
            else:
                visited_mols.add(mol)
        if not cyclic:
            filtered.append(p)
    return filtered

pathstat_f = open("pathfind-stats.txt", "w")

def find_gapless_fix(gapr, gapfix_state, params):#reachable_atoms, ag, scores, amp, atom_scores, path_reject_threshold, cdir):
    Q = priorityDictionary()
    Q[gapfix_state] = gapfix_state.cost

    visited = set()
    visited.add(hash(gapfix_state))
    visited_states = 0

    # block product atoms to prevent reverse of gapr to be used
    #params.forbidden_atoms = get_mol_atoms(amp.reactions[gapr].products, amp)
    params.forbidden_atoms = set()

    start_time = os.times()[0]
    complete_found = False

    results = []
    for state in Q:

        if len(results) >= MAX_RESULT_PATHWAYS:
            message(VERBOSE_SOME, "Max number of results reached")
            break

        if len(Q) > MAX_QUEUE_SIZE:
            if TERMINATE_ON_MAX_QUEUE:
                message(VERBOSE_SOME, "Max queue size exceeded")
                break
            else:
                continue

        now = os.times()[0]
        if now - start_time > MAX_SEARCH_TIME:
            message(VERBOSE_SOME, "Max search time exceeded")
            break

        message(VERBOSE_SOME, "%s (iter %d, state %d, %d in queue, %d results): %d r (%d paths), %d unsat, %d reached, %d used, %s cost" % (gapr, params.iteration, visited_states, len(Q), len(results), len(state.R), len(state.atom_paths), len(state.unsat), len(state.reached), len(state.used_atoms), state.cost))

        #if DRAW_PARTIAL_PATHWAYS:
        #    visualize_path(target_reactions, target_molecules, "partial%d" % (visited_states), state, sources, scores, amp, reachable_atoms, mol_names, atom_scores, cdir)

        visited_states += 1

        if visited_states > MAX_VISITED_STATES:
            message(VERBOSE_SOME, "Max number of partial pathways exceeded")
            break

        reaction_sets = {}

        #message(VERBOSE_ALL, "    KSP: k=%d: " % (NUM_PATHS)) 

        new_paths = find_atom_paths(NUM_PATHS, state.reached, state.unsat, params)

        if len(new_paths) == 0:
            #message(VERBOSE_ALL, "    KSP: unable to find paths")
            continue

        upaths = len(new_paths)
        new_paths = filter_paths(new_paths, params.amp)

        #message(VERBOSE_ALL, "%d/%d paths" % (upaths, len(new_paths)))

        reaction_sets = collect_reaction_sets(new_paths, state.unsat, params)

        #message(VERBOSE_ALL, "%d reaction sets" % (len(reaction_sets)))

        pathstat_f.write("%d\t%d\t%d\t%d\t%d\n" % (len(state.reached), len(state.unsat), upaths, len(new_paths), len(reaction_sets)))

        if len(reaction_sets) == 0:
            #message(VERBOSE_ALL, "    KSP: unable to find reactions")
            continue

        num_cyclic = 0
        num_already_done = 0
        num_too_costly = 0
        num_partial = 0
        num_complete = 0

        for adds in reaction_sets:
            new_state = state.copy()
            #cyclic = new_state.add_path(adds, amp, source_atoms, scores)
            cyclic = new_state.add_path(adds, params.amp, params.scores)

            if REJECT_CYCLIC and cyclic:
                tt = os.times()[0]
                #visualize_path([gapr], [], "%s-path" % (tt), new_state, set(), scores, amp, reachable_atoms, {}, atom_scores, odir, cdir)
                num_cyclic += 1
                continue

            # if COMPUTE_YIELD:
            #     new_state.tgt_yield = compute_yield(new_state, target_molecules, source_atoms, amp)
            #     if REJECT_NONPOSITIVE_YIELD and new_state.tgt_yield == 0:
            #         if verbose:
            #             print "        REJECT - NO TARGET YIELD"
            #         continue

            h = hash(new_state)
            if h in visited:
                num_already_done += 1
                continue
            visited.add(h)

            is_complete = new_state.complete()
            if is_complete:
                complete_found = True

            if new_state.total_log_prob() > params.path_reject_threshold * max_cost():
                num_too_costly += 1
                continue

            if not is_complete:
                est_cost = new_state.cost_estimate_indep(params.atom_scores)
                #if verbose:
                #    print "        Queue partial pathway", len(new_state.R), new_state.cost, est_cost
                Q[new_state] = est_cost
                num_partial += 1
            else:
                results.append(new_state)
                num_complete += 1
                #if verbose:
                #    print "        Complete pathway found", len(new_state.R), new_state.cost
                utime, stime, cutime, cstime, ertime = os.times()
                #outf.write("%s\t%s\t%d\t%s\n" % (gapr, utime, len(results), new_state.summarize()))

                if DRAW_COMPLETE_PATHWAYS:
                    #visualize_path([gapr], [], "%s-complete-%d" % (gapr, len(results)), new_state, sources, scores, amp, reachable_atoms, mol_names, atom_scores, odir, cdir)
                    pass
                if len(results) >= MAX_RESULT_PATHWAYS:
                    break
        #message(VERBOSE_ALL, "%d complete, %d partials, %d cyclic, %d already done, %d too costly" % (num_complete, num_partial, num_cyclic, num_already_done, num_too_costly))
    return results, complete_found
    
NO_FILL_FOUND = "NoGapFixFound"
ALREADY_IN_NET = "AlreadyInNetwork"
NOT_GAPPED = "NotGapped"
FILL_TOO_COSTLY = "FillTooCostly"
GAP_FILLED = "GapFilled"

def gapless_fix(gapr, params):
#, amp, ag, sources, ubiqs, scores, 
#                reachable_atoms, 
#                mol_names, atom_scores, path_reject_threshold,
#                odir, initial_state, source_atoms, cdir):

    current_state = params.current_state.copy()

    #print "*** Round %d: %d/%d: %d r, %d reached. Next: %s (cost %.2f) ***" % (iteration, gapix + 1, num_total, len(current_state.R), len(current_state.reached), gapr, rcost(scores[gapr].pscore))

    if gapr in current_state.R:
        message(VERBOSE_ALL, "%s already in network" % (gapr))
        return ALREADY_IN_NET, current_state

    # add next reaction in sequence to pathway, fix gaps that may result
    gapfix_state = SearchState()
    gapfix_state.reached = set(current_state.reached)

    re = params.amp.reactions[gapr]
    # ignore ubiquituous metabolites as reaction substrates
    sub_atoms = get_mol_atoms(re.substrates, params.amp)
    for u in sub_atoms:
        #print u, amp.get_mol_and_atom(u), u in ag.E, u in ag.Erev
        if u not in params.source_atoms and u in params.reachable_atoms and u in params.ag.E:
            gapfix_state.add_target(u)
    gapfix_state.add_reaction(gapr, params.amp, params.scores)

    if gapfix_state.complete():
        current_state.add_path(gapfix_state.R, params.amp, params.scores)
        message(VERBOSE_SOME, "%s not gapped" % (gapr))
        return NOT_GAPPED, current_state

    # search for a good gapfilling pathway for this reaction
    results, complete_found = find_gapless_fix(gapr, gapfix_state, params)
    results.sort(lambda x, y: cmp(x.cost, y.cost))
    utime, stime, cutime, cstime, ertime = os.times()

    if DRAW_CANDIDATE_PATHWAYS:
        for resix, res in enumerate(results):
            #visualize_path([gapr], [], "_%s-complete-%d" % (gapr, resix), res, sources, scores, amp, reachable_atoms, mol_names, atom_scores, odir, cdir)
            pass

    if len(results) > 0:
        best_state = results[0]
        best_state.update_cost(params.scores)
        gname = ""
        if best_state.total_log_prob() > params.path_reject_threshold * max_cost():
            message(VERBOSE_SOME, "%s: Best fill too expensive" % (gapr))
            gname = "rejected"
            #current_state.add_reaction(gapr, params.amp, params.scores) 
            return FILL_TOO_COSTLY, current_state
        else:
            message(VERBOSE_SOME, "%s: Gap fill found" % (gapr))
            gname = "accepted"
            current_state.add_path(best_state.R, params.amp, params.scores) 
            return GAP_FILLED, current_state

    else: # no results
        if complete_found:
            # a complete pathway found but too costly and not in results
            message(VERBOSE_SOME, "%s: Best gap fill too expensive" % (gapr))
            #current_state.add_reaction(gapr, params.amp, params.scores)
            return FILL_TOO_COSTLY, current_state
        else:
            message(VERBOSE_SOME, "%s: No gap fill found" % (gapr))
            #current_state.add_reaction(gapr, params.amp, params.scores) 
            return NO_FILL_FOUND, current_state

def choose_next_reaction(current_state, params):
    #print "choose_next_reaction: compute_atom_score_heuristic"
    reachable_atoms, atom_scores, preceding_reactions = compute_atom_score_heuristic(source_atoms = current_state.reached, params = params)
	
    current_base = set(map(lambda x: params.amp.reactions[x].basename(), current_state.R))
    candidates = []

    #print "current_base:", list(current_base)

    base_reactions = set(map(lambda x: params.amp.reactions[x].basename(), params.amp.reactions))

    if params.target_reactions != None:
        base_reactions.intersection_update(params.target_reactions)

    #print "Failed reactions:", params.failed_reactions
    #print "Current base:", current_base

    while len(candidates) == 0:
        
        for rbase in base_reactions:	   
            if rbase in params.failed_reactions or rbase in current_base:
                continue
            mapname = "%s_1" % (rbase)      # use reaction with atom map if available
            if mapname not in params.scores:
                mapname = "%s_0" % (rbase)  # use unmapped reaction instead
            if params.scores[mapname].pscore > params.threshold:
 		#print "candidate:", mapname	
                candidates.append(mapname)

        #print "Candidates:", candidates

        if len(candidates) == 0:
            params.iteration += 1
            #print "Iteration now", params.iteration
            if params.failed_reactions == params.prev_failed_reactions:
                print "Convergence"
                return None
            params.prev_failed_reactions = params.failed_reactions
            params.failed_reactions = set()

    max_re_dists = estimate_reaction_costs(atom_scores, candidates, params.amp)
    rord = max_re_dists.keys()
    rord.sort(lambda x, y: cmp(max_re_dists[x], max_re_dists[y]))

    #for r in rord:
    #    print r, max_re_dists[r]

    #print "Choosing", rord[0]

    return rord[0]

#def find_gapless_network(iteration, rord, amp, ag, sources, ubiqs, scores, 
#                         reachable_atoms, 
#                         mol_names, atom_scores, path_reject_threshold,
#                         odir, cdir):
def find_gapless_network(params):
    """
Reconstruct a gapless metabolic network by adding a single reaction at time and filling resulting gaps. 

Reaction addition order is determined by reaction distance from reachable metabolites.

When unable to add reaction, add it to failed list. After no more reactions can be added, consider all reactions on the failed list again. Terminate when process converges. Do not add failed reactions to the network.
"""

    outf = open("%s/%s" % (params.odir, common.FILLS_FILE), "w")
    outf.write("#Output of %s on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))
    outf.write("#TargetR EC Time Iteration Result UTime ")
    outf.write("%s\n" % (SearchState.summary_header()))

    params.source_atoms = get_mol_atoms(params.sources, params.amp)
    params.current_state = SearchState()
    params.current_state.reached = set(params.source_atoms)

    #stuff = params.amp.reactions.keys()[1:20]
    #for r in stuff:
    #    params.target_reactions.add(params.amp.reactions[r].basename())
    #    print "added", r
    
    decisions = {GAP_FILLED : [], 
                 NO_FILL_FOUND : [],
                 ALREADY_IN_NET : [],
                 FILL_TOO_COSTLY : [],
                 NOT_GAPPED : []}

    params.iteration = 1
    params.failed_reactions = set()
    params.prev_failed_reactions = None
    
    #for gapix, gapr in enumerate(params.rord):
    gapix = 0
    while 1:
        gapr = choose_next_reaction(params.current_state, params)
        if gapr == None:
            # No more reactions to add
            break
        gapix += 1

        fwd = "%s_1" % (gapr)
        rev = "%s_1_rev" % (gapr)
        
        if fwd not in params.amp.reactions:
            # only an unmapped reaction available, lets use that
            # TODO: how to clean up the choice of reaction
            fwd = "%s_0" % (gapr)
            rev = "%s_0_rev" % (gapr)

        # delete atom graph edges associated with the reverse direction
	bound = searchBound(gapr)

	
        params.ag.delete_reaction_edges(params.amp.reactions[rev])
        decision_fwd, state_fwd = gapless_fix(fwd, params) #amp, ag, sources, ubiqs, scores, reachable_atoms, mol_names, atom_scores, path_reject_threshold, odir, current_state, source_atoms, cdir)
        params.ag.add_reaction_edges(params.amp.reactions[rev], rcost(params.scores[rev].pscore))
	
        params.ag.delete_reaction_edges(params.amp.reactions[fwd])
        decision_rev, state_rev = gapless_fix(rev, params) #, amp, ag, sources, ubiqs, scores, reachable_atoms, mol_names, atom_scores, path_reject_threshold, odir, current_state, source_atoms, cdir)
        params.ag.add_reaction_edges(params.amp.reactions[fwd], rcost(params.scores[fwd].pscore))

        ACCEPTS = [ALREADY_IN_NET, NOT_GAPPED, GAP_FILLED]

        #print decision_fwd, state_fwd.cost, decision_rev, state_rev.cost
	
        if decision_rev not in ACCEPTS or (decision_fwd in ACCEPTS and state_fwd.cost <= state_rev.cost):
            #message(VERBOSE_ALL, "    Forward dir selected %s %s" % (state_fwd.cost, state_rev.cost))
            new_current_state = state_fwd
            decision = decision_fwd
        else:
            #message(VERBOSE_ALL, "    Reverse dir selected %s %s" % (state_fwd.cost, state_rev.cost))
            new_current_state = state_rev
            decision = decision_rev

        diff_state = new_current_state.copy()
        diff_state.R.difference_update(params.current_state.R)
        diff_state.update_cost(params.scores)
        utime, stime, cutime, cstime, ertime = os.times()
        outf.write("%s\t%s\t%s\t%d\t%s\t%s\t%s\n" % (gapr, ",".join(params.scores[fwd].ec), datetime.datetime.now(), params.iteration, decision, utime, diff_state.summarize()))
        outf.flush()
        #visualize_path([gapr], [], "%s-%s" % (gapr, decision), diff_state, params) #sources, scores, amp, reachable_atoms, mol_names, atom_scores, odir, cdir)

        decisions[decision].append(gapr)

        if decision == NO_FILL_FOUND or decision == FILL_TOO_COSTLY:
            params.failed_reactions.add(gapr)

        params.current_state = new_current_state

    if ADD_FAILED_HIGHSCORE_REACTIONS:
        for r in params.failed_reactions:
            rn = "%s_1" % (r)
            if rn not in params.amp.reactions:
                rn = "%s_0" % (r)
            message(VERBOSE_ALL, "Adding high-scoring failed reaction %s" % (rn))
            params.current_state.add_reaction(rn, params.amp, params.scores)

    oo = open("%s/%s" % (params.odir, common.SPONTANEOUS_FILE), "w")

    if ADD_SPONTANEOUS_REACTIONS:
        for br in params.amp.spontaneous_reactions:
            for r in params.amp.base_reactions[br]:
                #rn = "%s_0" % (r)
                message(VERBOSE_ALL, "Adding spontaneous reaction %s" % (r))
                params.current_state.add_reaction(r, params.amp, params.scores)
                oo.write("""%s	%s	?	?	?	Spontaneous	0.02	0.0	0	0.0	0.0	0.0	None	None	0.0	None	None\n""" % (br, r))

    return params.current_state, decisions

def calculate_max_mol_distances(atom_scores, amp):
    max_dists = {}
    dists = {}
    for u in atom_scores:
        mol, atom = amp.get_mol_and_atom(u)
        if mol not in dists:
            dists[mol] = []
        dists[mol].append(atom_scores[u])
        if mol not in max_dists or atom_scores[u] > max_dists[mol]:
            max_dists[mol] = atom_scores[u]
 
    #for mol in dists:
        #dists[mol].sort()
        #print mol, ",".join(map(str, dists[mol]))

    return max_dists    

def gapless(amp, net, sources, ubiqs, scores, mol_names, threshold, path_reject_threshold, odir, ecscores, r2ec, cdir, target_reactions = None):

    sys.setrecursionlimit(100000) # for Tarjan's algorithm in scc.py

    params = Parameters()
    params.amp = amp
    params.net = net
    params.sources = sources
    params.ubiqs = ubiqs
    params.scores = scores
    params.mol_names = mol_names
    params.threshold = threshold
    params.path_reject_threshold = path_reject_threshold
    params.odir = odir
    params.ecscores = ecscores
    params.r2ec = r2ec
    params.cdir = cdir

    message(VERBOSE_SOME, "Building atom graph...")
    params.ag = network.AtomGraph(amp.reactions, lambda x: rcost(scores[x].pscore))

    message(VERBOSE_SOME, "Scaling atom graph edge weights...")
    params.ag.scale_edge_weights(amp)

    message(VERBOSE_SOME, "Computing reach heuristic...")
    params.reachable_atoms, params.atom_scores, params.preceding_reactions = compute_atom_score_heuristic(sources = params.sources, params =  params)
    message(VERBOSE_SOME, "%d reachable atoms" % (len(params.reachable_atoms)))

    params.reachable_mols = set()
    for u in params.reachable_atoms:
        params.reachable_mols.add(amp.get_mol_and_atom(u)[0])

    params.target_reactions = target_reactions

    #message(VERBOSE_SOME, "Computing max mol distances...")
    #max_mol_dists = calculate_max_mol_distances(atom_scores, amp)

    #print "Computing max reaction distances..."

    write_parameters(params)

    fixed, decisions = find_gapless_network(params)

    netf = open("%s/%s" % (odir, common.NETWORK_REACTION_FILE), "w")
    netf.write("#Generated by %s in %s\n" % (sys.argv[0], datetime.datetime.now()))
    res = list(fixed.R)
    res.sort()
    netf.write("#Reaction ECs PScore Cost Level %s\n" % (scores[res[0]].header()))
    for r in res:
        base = amp.reactions[r].basename()
        score = scores[r]
        if base in r2ec:
            ecs = r2ec[base]
        else:
            ecs = ["?"]
        
        if score.pscore > threshold:
            gaps = "High"
        else:
            gaps = "Low"

        netf.write("%s\t%s\t%f\t%f\t%s\t%s\n" % (r, ",".join(ecs), score.pscore, rcost(score.pscore), gaps, score))

        for ec in ecs:
            if ec not in ecscores or ecscores[ec].pscore < score.pscore:
                ecscores[ec] = score

    netf.close()

    ecf = open("%s/%s" % (odir, common.NETWORK_EC_FILE), "w")
    keys = ecscores.keys()
    keys.sort()
    ecf.write("#Generated by %s in %s\n" % (sys.argv[0], datetime.datetime.now()))
    ecf.write("#EC PScore Cost Level %s\n" % (ecscores[keys[0]].header()))
    for ec in keys:
        score = ecscores[ec]
        if score.pscore > threshold:
            ecf.write("%s\t%s\t%s\t%s\t%s\n" % (ec, score.pscore, rcost(score.pscore), "High", score))
    ecf.close()

    #message(VERBOSE_SOME, "Generating SBML...")
    #network2sbml.main(odir, common.FILE_STOICHIOMETRY, common.FILE_METABOLITE_NAMES, "%s/%s" % (odir, common.FILE_SBML_OUTPUT))

    ## New reaction bag doesn't contain kegg-pathway mappings 
    ## Skipping this tep entirely
    #message(VERBOSE_SOME, "Mapping pathways to KEGG...")
    #map_result_in_kegg.main(odir, odir)

    message(VERBOSE_SOME, "\o/ All done \o/")

def main(bfile, cdir, kdir, adir, sourcesfn, ubiqfn, scoresfn, threshold, path_reject_threshold, odir, cnamefile, ecfile, target_reactions = None):   
    readBoundaries(bfile)
    print bfile
    try:
        os.mkdir(odir)
    except:
        pass

    ec2r, r2ec = common.read_ec_list(open(ecfile))

    precalc_colors()
    message(VERBOSE_SOME, "Loading molecule names... ")
    #  dictionary where keys are mol ids (C00001) and
    # items are names from second column of kegg-compounds file
    mol_names = common.parse_molecule_names(cnamefile)

    message(VERBOSE_SOME, "Loaded %d molecules." % (len(mol_names)))
    message(VERBOSE_SOME, "Loading atom maps... ")
    amp = atommap.AtomMappingParser()
    amp.readBoundaries(bfile)
    #accepted_atom_types = set(["C"])
    accepted_atom_types = set(DEFAULT_ATOM_TYPES)
    skip_atom_types = None
    #skip_atom_types = set(["O", "H"])
    #skip_atom_types = set(["O", "H", "N", "P", "R", "S", "Mn"])
    amp.parse_reaction_dir(adir, accepted_atom_types, skip_atom_types)
    message(VERBOSE_SOME, "Loaded %d mapped reactions." % (len(amp.reactions)))
    message(VERBOSE_SOME, "Loading KEGG LIGAND...")
    amp.parse_ligand_reaction(kdir)
    message(VERBOSE_SOME, "Loaded %d reactions in total." % (len(amp.reactions)))
    message(VERBOSE_SOME, "Loading source metabolites... ")
    sources = common.read_set(open(sourcesfn))
    #print sources	
    message(VERBOSE_SOME, "Loaded %d sources." % (len(sources)))
    message(VERBOSE_SOME, "Loading ubiquituous metabolites... ")
    ubiqs = common.read_set(open(ubiqfn))
    message(VERBOSE_SOME, "Loaded %d metabolites." % (len(ubiqs)))
    message(VERBOSE_SOME, "Loading reaction scores... ")
    sys.stdout.flush()
    scores, ecscores = read_reaction_scores(scoresfn, amp, r2ec, threshold)
    message(VERBOSE_SOME, "Loaded %d scored reactions." % (len(scores)))
    #scores = generate_scores(amp.reactions)
    message(VERBOSE_SOME, "Setting up universal metabolic network... ")
    sys.stdout.flush()
    net = network.MetabolicNetwork(amp.reactions)

    gapless(amp, net, sources, ubiqs, scores, mol_names, threshold, path_reject_threshold, odir, ecscores, r2ec, cdir, target_reactions)

if __name__ == "__main__":
    bfile = sys.argv[1]	     # bounds file	
    kdir = sys.argv[2]       # KEGG LIGAND dir
    adir = sys.argv[3]       # dir with atom maps
    sourcesfn = sys.argv[4]  # list of source metabolites
    ubiqfn = sys.argv[5]     # list of ubiquituous metabolites 
    scoresfn = sys.argv[6]   # reaction scores
    threshold = float(sys.argv[7])
    path_reject_threshold = float(sys.argv[8])
    odir = sys.argv[9]       # output dir
    cnamefile = sys.argv[10]     # list of ubiquituous metabolites 
    ecfile = sys.argv[11]
    if len(sys.argv) > 12:
        target_reactions = set(sys.argv[12].split(","))
    else:
        target_reactions = None

    cdir = "."
    main(bfile, cdir, kdir, adir, sourcesfn, ubiqfn, scoresfn, threshold, path_reject_threshold, odir, cnamefile, ecfile, target_reactions)
