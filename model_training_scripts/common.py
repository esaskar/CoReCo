"""Global parameters and convience functions.

"""

import os, math, datetime

import tree

# Global parameters

DATE_FORMAT = "%Y.%m.%d %H:%M:%S"

# Sequence analysis parameters

BLAST_PVALUE_CUTOFF = 1e-200
EXP_APPROX = 1e-5

# Global input files

AUX_DIR = "../../data/Kegg/aux"
#FILE_EC_MAP = "ec-list.txt"
#FILE_EC_MAP = "%s/ec-list-augmented.txt" % (AUX_DIR) # contains manual fixes to KEGG mappings; see README
FILE_EC_MAP = "%s/ec-list.txt" % (AUX_DIR) # contains manual fixes to KEGG mappings; see README
FILE_METABOLITE_NAMES = "%s/kegg-compounds" % (AUX_DIR) 
FILE_STOICHIOMETRY = "kegg/balances.eqn"
FILE_SBML_OUTPUT = "network.sbml"


# Reconstruction project files

PATHWAYS_DIR = "pathways"
FILLS_FILE = "fills"
NETWORK_EC_FILE = "network.ecs"
NETWORK_ECLIST_FILE = "network.eclist"
NETWORK_REACTION_FILE = "network.reactions"
NETWORK_FILE = "network"

REACTION_SCORE_DIR = "rscores"
RECO_RESULT_DIR = "reco"
FILE_EC_GRAPH = "network.ecgraph"

FILE_ECS_TP = "comparison-ec.tps"
FILE_ECS_FP = "comparison-ec.fps"
FILE_ECS_FN = "comparison-ec.fns"

PARAMETER_FILE = "parameters"

EVIDENCE_FILE = "evidence"
REACTION_SCORE_FILE = "reaction-scores"
FULL_REACTION_SCORE_FILE = "reaction-scores.full"
ROC_CURVE_DIR = "rocs"
CPD_DIR = "cpds"
ORGANISM_LIST_FILE = "org_list.backup"
SOURCE_METABOLITE_FILE = "sources"
REACTION_SCORE_PLOT_DIR = "score-plots"
EC_VISUALS_DIR = "ec-visuals"
TREE_CPD_FILE = "tree.cpds"

def floatNone(x):
    if x == "None":
        return 0
    else:
        return float(x)

class Score(object):
    def __init__(self, s = None):        
        self.pscore = None
        self.npscore = None
        self.logratio = None
        self.btscore = None
        self.gtscore = None
        self.bscore = None
        self.bseq1 = self.bseq2 = None
        self.gscore = None
        self.gseq1 = self.gseq2 = None
        if s != None:
            vals = s.strip().split("\t")
            if len(vals) == 12:
                # fixes bad format produced by network2matrix.py
                vals = vals[1:]
            if len(vals) == 5:
                self.pscore, self.npscore, self.logratio, self.btscore, self.gtscore = map(float, vals)
            elif len(vals) == 11:
                self.pscore = floatNone(vals[0])
                self.npscore = floatNone(vals[1])
                self.logratio = floatNone(vals[2])
                self.btscore = floatNone(vals[3])
                self.gtscore = floatNone(vals[4])
                if vals[5] != "?":
                    self.bscore = float(vals[5])
                else:
                    self.bscore = 0.0
                self.bseq1 = vals[6]
                self.bseq2 = vals[7]
                if vals[8] != "?":
                    self.gscore = float(vals[8])
                else:
                    self.gscore = 0.0
                self.gseq1 = vals[9]
                self.gseq2 = vals[10]
            else:
                raise Exception("Parse error: expecting 5 or 11 values in score file: %s" % (s))

    @classmethod
    def strheader(cls, full = True):
        if full:
            return "FullP NaiveP PLogRatio BlastTreeP GTGTreeP BlastScore BSeq1 BSeq2 GTGScore GSeq1 GSeq2"
        else:
            return "FullP NaiveP PLogRatio BlastTreeP GTGTreeP"

    def header(self):
        if self.bscore == None:
            return "FullP NaiveP PLogRatio BlastTreeP GTGTreeP"
        else:
            return "FullP NaiveP PLogRatio BlastTreeP GTGTreeP BlastScore BSeq1 BSeq2 GTGScore GSeq1 GSeq2"

    def __str__(self):
        if self.bscore == None:
            return "%s\t%s\t%s\t%s\t%s" % (self.pscore, self.npscore, self.logratio, self.btscore, self.gtscore)
        else:
            return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.pscore, self.npscore, self.logratio, self.btscore, self.gtscore, self.bscore, self.bseq1, self.bseq2, self.gscore, self.gseq1, self.gseq2)

class ScoredReaction:
    def __init__(self, rid, ec, p, cost, level, score):
        self.rid = rid
        #self.rmap = rmap
        self.ec = ec
        self.p = p
        self.cost = cost
        self.level = level
        self.score = score
    @classmethod
    def strheader(cls):
        return "Reaction EC PScore Cost Level %s" % (Score.strheader())
    def __repr__(self):
        return "%s:%s" % (self.rid, self.p)
    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s" % (self.rid, self.ec, self.p, self.cost, self.level, self.score)

class ReactionWrapper:
    def __init__(self, reactions):
        self.reactions = reactions

class FillReaction:
    def __init__(self, values):
        #TargetR EC Time Iteration Result UTime NumReactions TargetYield TotalCost AvgCost TotalLogP TotalP AvgLogP NumPaths Reactions
        self.rid, self.ec, self.recotime, self.iteration, self.decision, self.utime, self.nreactions, self.tyield, self.tcost, self.avgcost, self.tlp, self.totalp, self.alp, self.npaths = values[0:14]
        if len(values) == 15:
            self.reactions = values[14].split(",")
        else:
            self.reactions = [self.rid]
        self.ec = self.ec.split(",")
        #self.recotime = float(self.recotime)
        self.utime = float(self.utime)
        self.nreactions = int(self.nreactions)
        self.tcost = float(self.tcost)
        self.avgcost = float(self.avgcost)
        self.tlp = float(self.tlp)
        self.totalp = float(self.totalp)
        self.alp = float(self.alp)
        self.npaths = int(self.npaths)

def log(msg):
    print "%s: %s" % (datetime.datetime.now().strftime(DATE_FORMAT), msg)

def read_lines(f):
    for s in f:
        if s.startswith("#"):
            continue

        yield s.strip().split("\t")

def read_set(f):
    a = set()
    for s in read_lines(f):
        a.add(s[0])
    return a

def write_set(f, A):
    S = list(A)
    S.sort()
    for a in S:
        f.write("%s\n" % (a))

def read_dict(f, vsplit = None):
    D = {}
    for k, v in read_lines(f):
        if vsplit:
            v = v.split(vsplit)
        D[k] = v
    return D

def read_source_metabolites(ddir):
    f = open("%s/%s" % (ddir, SOURCE_METABOLITE_FILE))
    return read_set(f)

def read_organisms(ddir):
    orgnames = {}
    longtoshort = {}
    f = open("%s/%s" % (ddir, ORGANISM_LIST_FILE))
    for s in f:
        short, lname = s.strip().split()
        orgnames[short] = lname
        longtoshort[lname] = short
    return longtoshort, orgnames

def score_transform(x):
    return math.log(x) / math.log(2)

def read_model(f, filter_partial_ecs = False):
    model = set()
    for s in f:
        if s.startswith("#"):
            continue
        ec = s.strip()
        if filter_partial_ecs and "-" in ec:
            continue
        model.add(ec)
    return model

def read_models(mdir, filter_partial_ecs = False):
    models = {}
    fns = os.listdir(mdir)
    ecs = set()
    for fn in fns:
        models[fn] = read_model(open("%s/%s" % (mdir, fn)), filter_partial_ecs)
        ecs.update(models[fn])
    return models, ecs

def read_scores(ddir, augmented_format = False):
    if augmented_format:
        fn = FULL_REACTION_SCORE_FILE
    else:
        fn = REACTION_SCORE_FILE
    f = open("%s/%s" % (ddir, fn))
    scores = {}
    ecs = set()
    for s in f:
        if s.startswith("#"):
            continue
        ec, species, scorestr = s.strip().split("\t", 2)
        score = Score(scorestr)
        ecs.add(ec)
        if species not in scores:
            scores[species] = {}
        scores[species][ec] = score
    return scores, ecs

def read_tree(fn):
    return tree.newickToTree(open(fn).read())

def read_tree_adj_list(ddir):
    fn = "%s/tree.full" % (ddir)
    f = open(fn)
    parent, children = {}, {}
    for s in f:
        if s.startswith("#"):
            continue
        u, v = s.strip().split("\t")
        if u not in children:
            children[u] = []
            parent[u] = None
        if v not in children:
            children[v] = []
            parent[v] = None
        children[u].append(v)
        parent[v] = u
    return parent, children


def read_fills(ddir):
    f = open("%s/%s" % (ddir, FILLS_FILE))
    results = {}
    for line in read_lines(f):
        results[line[0]] = FillReaction(line)
    return results
        
def read_ec_list(f):
    ec2r = {}
    r2ec = {}
    for s in f:
        if s.startswith("#"):
            continue
        ec, rs = s.strip().split("\t")
        rs = rs.split(",")
        ec2r[ec] = set(rs)
        for r in rs:
            if r not in r2ec:
                r2ec[r] = set()
            r2ec[r].add(ec)
    return ec2r, r2ec

def read_reconstruction(f):
    reco = {}
    for vals in read_lines(f):
        rid, ec, fullp, cost, level = vals[:5]
        if len(vals) > 5:
            score = Score("\t".join(vals[5:]))
        else:
            score = None
        reco[rid] = ScoredReaction(rid, ec, fullp, cost, level, score)
    return reco

def parse_reactants(s):
    mols = s.strip().split(" + ")
    stuff = set()
    coeff = {}
    for mol in mols:
        vals = mol.split(" ")
        if len(vals) == 2:
            qty, mol = vals
            qty = float(qty)
        else:
            qty = 1
        stuff.add(mol)
        coeff[mol] = qty
    return stuff, coeff

def read_stoichiometry(f):
    from atommap import Reaction
    reactions = {}
    for s in f:
        if s.startswith("#"):
            continue
        rid, balanced, origstatus, newstatus, eqn = s.strip().split("\t")
        if balanced != "True":
            #print "Skip %s" % (rid)
            continue
        lhs, rhs = eqn.split(" <=> ")
        subs, sub_coeff = parse_reactants(lhs)
        pros, pro_coeff = parse_reactants(rhs)
        r = Reaction(rid, subs, pros, sub_coeff, pro_coeff, {}, {})
        reactions[rid] = r
        rev = r.reverse()
        reactions[rev.id] = rev

    return ReactionWrapper(reactions)

def read_metabolite_names(f):
    names = {}
    for vals in read_lines(f):
        names[vals[0]] = vals[1]
    return names
        
def log(msg):
    print "%s: %s" % (datetime.datetime.now().strftime(DATE_FORMAT), msg)

def blast_evalue_to_pvalue(ev):
    if float(ev) > EXP_APPROX:
        pv = 1 - math.e ** (-ev)
    else:
        pv = ev
    return pv

def blast_joint_score(fwdp, revp):
    if fwdp < BLAST_PVALUE_CUTOFF:
        fwdp2 = BLAST_PVALUE_CUTOFF
    else:
        fwdp2 = fwdp
    if revp < BLAST_PVALUE_CUTOFF:
        revp2 = BLAST_PVALUE_CUTOFF
    else:
        revp2 = revp
    return -math.log(fwdp2 + revp2 - fwdp2 * revp2) / math.log(2)

def max_blast_joint_score():
    return blast_joint_score(BLAST_PVALUE_CUTOFF, BLAST_PVALUE_CUTOFF)

if __name__ == "__main__":
    pass
