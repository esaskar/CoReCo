#!/usr/bin/env python

import sys, datetime
from optparse import OptionParser

import atommap, common

import glpk

PO_CONSTRAINT = 1

ALL_UNIDIRECTIONAL = 0

# kJ/mol
GIBBS_THRESHOLD = 5.0

FBA_EPS = 1e-6

PRINT_FBA = 0
PRINT_MATRICES = 0

STATUS_OPTIMAL = 1

DEADEND_SOURCES = 0
DEADEND_INTAKE = 1.0

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

# Maximize
#  obj: + 0.6 x1 + 0.5 x2

# Subject To
#  c1: + x1 + 2 x2 <= 1
#  c2: + 3 x1 + x2 <= 2

# Bounds
#  x1 free
#  x2 free

# End

def lpname(n):
    return n.replace("-", "_")

def write_lp(f, N, objix, intakes, bounds, ix2m, ix2r):
    f.write("Maximize\n")
    s = " obj: + 1.0 %s" % (lpname(ix2r[objix]))
    f.write("%s\n" % (s))

    f.write("Subject To\n")
    for mix in ix2m:
        s = " %s:" % (lpname(ix2m[mix]))
        for rix in N:
            if mix in N[rix]:
                v = N[rix][mix]
                if v != 0:
                    s += " %+f %s" % (v, lpname(ix2r[rix]))
        s += " = 0"
        f.write("%s\n" % (s))

    f.write("Bounds\n")
    for rix in intakes:
        r = ix2r[rix]
        if intakes[rix] == 0:
            f.write(" %s = 0\n" % (lpname(r)))
        else:
            f.write(" %+e <= %s <= %+e\n" % (0, lpname(r), intakes[rix]))
    for rix in ix2r:
        r = ix2r[rix]
        if r in bounds:
            lb, ub = bounds[r]
            if lb == None:
                lbs = ""
            else:
                lbs = "%+e <=" % (lb)
            if ub == None:
                ubs = ""
            else:
                ubs = "<= %+e" % (ub)
            if lb == 0 and ub == 0:
                f.write(" %s = 0\n" % (lpname(r)))
            else:
                f.write(" %s %s %s\n" % (lbs, lpname(r), ubs))        

    f.write("End\n")
    f.close()


def compute_yield_multiprocessing(q, network, objective, source_mols, amp, mol_names, options):
    ret = compute_yield(network, objective, source_mols, amp, mol_names, options)
    q.put(ret)

def compute_yield_process(model, objective, src, amp, mol_names, options):
    import multiprocessing
    q = multiprocessing.Queue()
    p = multiprocessing.Process(target=compute_yield_multiprocessing, args=(q, model, objective, src, amp, mol_names, options))
    p.start()
    p.join()
    return q.get()

def compute_yield(network, objective, source_mols, amp, mol_names, options):

    #network = set(list(network)[:1])
    #print network

    #objective = {}
    #objective["C14101"] = 1
    #source_mols = {"C00001" : 1.0, "C14120" : 1.0}

    if ALL_UNIDIRECTIONAL:
        for r in amp.reactions:
            if r in options.bounds:
                continue
            if r.endswith("_rev"):
                options.bounds[r] = (0, 0)
            else:
                options.bounds[r] = (0, None)

    all_mols = set(objective.keys())
    for r in network:
        re = amp.reactions[r]
        all_mols.update(re.substrates)
        all_mols.update(re.products)
        #print r, re.substrates, re.products

    POC = "P/O"
    if PO_CONSTRAINT:
#        print "#!!! P/O constraint hardcoded !!!"
        all_mols.add(POC)
        amp.reactions["R00081"].products.add(POC)
        amp.reactions["R00086"].products.add(POC)
        amp.reactions["R00081"].pro_coeff[POC] = 2
        amp.reactions["R00086"].pro_coeff[POC] = 1
        amp.reactions["R00081_rev"].substrates.add(POC)
        amp.reactions["R00086_rev"].substrates.add(POC)
        amp.reactions["R00081_rev"].sub_coeff[POC] = 2
        amp.reactions["R00086_rev"].sub_coeff[POC] = 1

#        print "#!!! Added reactions: !!!"
#        print "#", amp.reactions["R00081"].description()
#        print "#", amp.reactions["R00086"].description()
#        print "#", amp.reactions["R00081_rev"].description()
#        print "#", amp.reactions["R00086_rev"].description()

    #print source_mols
    #print objective

    #import random
    #random.seed(0)
    #objective[random.choice(list(all_mols))] = -1

    #print "RANDOM OBJECTIVE", objective

    #print "SOURCES", source_mols

    sources = {}
    for mol in source_mols:
        if mol in all_mols:
            sources[mol] = source_mols[mol]

    #rint "%d mols" % (len(all_mols))
    #rint "%d source mols" % (len(sources))

    #no_intake_mols = all_mols.difference(source_mols)

    #all_mols = set(balanced_mols).union(source_mols)
    #not_balanced = all_mols.difference(balanced_mols)

    N = {} # stoichiometric matrix

    OBJ_INDEX = 0   # objective function index
    N[OBJ_INDEX] = {}

    if PRINT_FBA:
        print "OBJECTIVE -> %d" % (OBJ_INDEX)

    r2ix = {}
    m2ix = {}

    r2ix["OBJECTIVE"] = OBJ_INDEX

    global next_mol_index
    next_mol_index = 0
    next_r_index = OBJ_INDEX + 1
    for r in network:
        N[next_r_index] = {}
        r2ix[r] = next_r_index
        if PRINT_FBA:
            print r, "->", next_r_index
        re = amp.reactions[r]

        for mol in re.substrates:
            mol_ix = get_index(mol, m2ix)
            if mol_ix not in N[next_r_index]:
                N[next_r_index][mol_ix] = 0
            N[next_r_index][mol_ix] += -re.sub_coeff[mol]
            #print r, re.sub_coeff
            assert(-re.sub_coeff[mol] <= 0)
        for mol in re.products:
            mol_ix = get_index(mol, m2ix)
            if mol_ix not in N[next_r_index]:
                N[next_r_index][mol_ix] = 0
            N[next_r_index][mol_ix] += re.pro_coeff[mol]
            assert(re.pro_coeff[mol] >= 0)
        next_r_index += 1

    # add metabolite sink reactions and set up objective vector
    for mol in all_mols:
        if PO_CONSTRAINT and mol == POC:
#            print "#!!! Skipping POC SINK - POC metabolite in steady-state !!!"
            continue

        mol_ix = get_index(mol, m2ix)
        if PRINT_FBA:
            print "SINK", mol, mol_ix, "->", next_r_index
        r2ix["SINK-%s" % (mol)] = next_r_index
        N[next_r_index] = {}
        N[next_r_index][mol_ix] = -1.0
        if mol in objective:
            N[OBJ_INDEX][mol_ix] = objective[mol]
            #print "OBJ", N[OBJ_INDEX][mol_ix]
        #if mol == "C00025":
        #    print OBJ_INDEX, mol_ix, N[OBJ_INDEX][mol_ix]
        #    print next_r_index, mol_ix, N[next_r_index][mol_ix]
        next_r_index += 1

    # add metabolite source reactions
    intakes = {}
    for mol in sources:
        if sources[mol] == 0:
            continue
        mol_ix = get_index(mol, m2ix)
        r2ix["SOURCE-%s" % (mol)] = next_r_index
        if PRINT_FBA:
            print "SOURCE", mol, mol_ix, "->", next_r_index
        N[next_r_index] = {}
        N[next_r_index][mol_ix] = 1.0
        #if mol == "C00025":
        #    print next_r_index, mol_ix, N[next_r_index][mol_ix], sources[mol]
        intakes[next_r_index] = sources[mol]
        next_r_index += 1

    #for m in N:
    #    for n in N[m]:
    #        print m, n, N[m][n]

    ix2r = rev_index(r2ix)
    ix2m = rev_index(m2ix)
    num_r = next_r_index
    num_m = next_mol_index

    # set up linear program for FBA

    # objective                       -> c
    # reaction rates 0 <= x <= x_max  -> G + h
    # steady-state                    -> A + b

    lp = glpk.LPX()
    lp.name = "fba"
    lp.obj.maximize = True
    #lp.obj.maximize = False
    lp.rows.add(num_m) # metabolites
    lp.cols.add(num_r) # reactions

    for i in range(num_m):
        lp.rows[i].bounds = 0.0, 0.0  # this does not work with GLPK!
        #lp.rows[i].bounds = -1e-9, 1e-9  # steady state
        #lp.rows[i].bounds = 0, 1e-9  # steady state
        lp.rows[i].name = ix2m[i]
        if PRINT_FBA:
            print "STEADY-STATE", i, ix2m[i]

    boundsstr = []
    for i in range(num_r):
        rid = ix2r[i]
        lp.cols[i].name = rid

        if rid in options.bounds:
            lb, ub = options.bounds[rid]
            #assert(not (lb == 0 and ub == 0))
            lp.cols[i].bounds = lb, ub
            boundsstr.append("%s %s %s" % (rid, lb, ub))
        else:
            lp.cols[i].bounds = 0, None   # fluxes >= 0

        if PRINT_FBA:
            print "reaction", lp.cols[i].name, lp.cols[i].bounds

    boundsstr.sort()
    #print "=== FBA bounds ==="
    #print "\n".join(boundsstr)
    #print "=================="

    for i in intakes:
        if intakes[i] == 0:
            lp.cols[i].bounds = 0, 1e-9
            lp.cols[i].bounds = 0, 0
        else:
            lp.cols[i].bounds = 0, intakes[i]
        if PRINT_FBA:
            print "INTAKE", i, intakes[i]

    #print "!!! C00282 output constrained !!!"
    #lp.cols[r2ix["SINK-C00282"]].bounds = 0, 1e-2

    objv = [0 for x in range(num_r)]
    objv[OBJ_INDEX] = 1.0

    if 0:
        print "#!!! maximizing R00086_rev !!!"
        objv[OBJ_INDEX] = 0.0
        objv[r2ix["R00086_rev"]] = 1.0
        lp.cols[r2ix["R00086"]].bounds = 0, 0
        lp.cols[r2ix["R00086_rev"]].bounds = 100, 100
        #objv[r2ix["R00081"]] = 1.0
        #lp.cols[r2ix["R00081"]].bounds = 100, 1000
        #lp.cols[r2ix["R00081_rev"]].bounds = 0, 0

#    print "!!! OBJECTIVE BOUNDED !!!"
#    lp.cols[OBJ_INDEX].bounds = 0, 100000.0

    #write_lp(open("lp", "w"), N, OBJ_INDEX, intakes, options.bounds, ix2m, ix2r)
    #sys.exit()

    lp.obj[:] = objv
    M = []

    for ri in N:
        for mi in N[ri]:
            M.append((mi, ri, N[ri][mi]))

    #print "Reactions:", num_r
    #print "Metabolites:", num_m

    lp.matrix = M

    o = open("fba-fluxes", "w")

    retval = lp.simplex()
    if lp.status != 'opt' and lp.status != "undef" and lp.status != "unbnd": 
        print lp.obj.value
        print lp.status
        print retval
        del lp
        return None, "no-solution"

    if lp.status != "opt":
        sys.stderr.write("#Warning: solution status: %s\n" % (lp.status))

    bal = 0.0
    tgt = "C00002"

    obj_comp_yields = {}
    for mol in objective:
        obj_comp_yields[mol] = 0.0

    allzero = True
    import math
    for c in lp.cols:
        if abs(c.primal) > 1e-7:
            allzero = False
            try:
                na = c.name.split("-")[1]
            except:
                na = ""
            if na in mol_names:
                na = mol_names[na]
            if c.name in amp.reactions:
                eqn = "%s <=> %s" % (map(mol_names.get, amp.reactions[c.name].substrates), map(mol_names.get, amp.reactions[c.name].products))
                eqn = "%s <=> %s" % (" + ".join(amp.reactions[c.name].substrates), " + ".join(amp.reactions[c.name].products))

                r = amp.reactions[c.name]
                for sub in r.sub_coeff:
                    if sub in obj_comp_yields:
                        obj_comp_yields[sub] -= r.sub_coeff[sub] * c.primal
                for pro in r.pro_coeff:
                    if pro in obj_comp_yields:
                        obj_comp_yields[pro] += r.pro_coeff[pro] * c.primal

            else:
                eqn = ""
            if options.fluxes and ("SINK" in c.name or "SOURCE" in c.name or options.verbose):
                if c.primal == c.bounds[1]:
                    limited = "Constrained"
                else:
                    limited = ""
                #print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (c.name, c.primal, na, eqn, c.bounds[0], c.bounds[1], limited)
                print "%s\t%s\t%s\t%s\t%s" % (c.name, c.primal, na, eqn, limited)

            if abs(c.primal) > 1e-10:
                o.write("%s\t%s\n" % (c.name, c.primal))
    #if allzero:
    #    print "Trivial solution"
    #exit()

    ov = lp.obj.value
    del lp

    #print "Balance of %s = %s" % (tgt, bal)
    for mol in obj_comp_yields:
        if mol in mol_names:
            na = mol_names[mol]
        else:
            na = "?"
        #print "OBJ-%s\t%s\t%s" % (mol, obj_comp_yields[mol], na)

    return ov, STATUS_OPTIMAL

def compute_yield_from_state(state, target_molecules, sources, amp, mol_names, options):
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

    objective = {}
    for mol in target_molecules:
        objective[mol] = 1

    return compute_yield(state.R, objective, source_mols, balanced_mols, amp, mol_names, options)

def test_objective_components(model, objective, sources, amp, of, mol_names, options):
    for mol in objective:
        obj = {mol : -1}
        res, status = compute_yield(model, obj, sources, amp, mol_names, options)
        of.write("%s\t%s\n" % (mol, res))
        of.flush()

def compare_models(model1, model2, objective, sources, amp, of, mol_names, options):
    of.write("#Mol MolName Model1 Model2 Diff Frac\n")
    for mol in objective:
        obj = {mol : -1}
        res1, status1 = compute_yield(model1, obj, sources, amp, mol_names, options)
        if model2 != None:
            res2, status2 = compute_yield(model2, obj, sources, amp, mol_names, options)
            diff = res2 - res1
            EPS = 1e-6
            frac = (res1 + EPS) / (res2 + EPS)
        else:
            res2 = float("nan")
            diff = float("nan")
            frac = float("nan")
        of.write("%s\t%s\t%f\t%f\t%f\t%f\n" % (mol, mol_names[mol], res1, res2, diff, frac))
        of.flush()

def is_number(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

def check_reactant_coefficients(model, amp):
    for r in model:
        re = amp.reactions[r]
        not_num = False
        for mol in re.substrates:
            if not is_number(re.sub_coeff[mol]):
                not_num = True
                break
        for mol in re.products:
            if not is_number(re.pro_coeff[mol]):
                not_num = True
                break
        if not_num:
            print r, "non-numeric coeff", re.sub_coeff, re.pro_coeff

import re
# "(2n+1)"
re_coeff = re.compile("(\(?[\d*nm]?[\+\-]?[\d*nm]?\)?)")
#"(12n)"
re_coeff2 = re.compile("(\d+)(n)")
# "C00031(n+1)"
re_coeff3 = re.compile("(\w\d{5})\(([mn\d*]?[\+\-]?[mn\d*]?)\)")

def convert_coeff(x):
    c = re_coeff.findall(str(x))[0]
    m = re_coeff2.match(c)  # e.g. "2n"
    if m != None:
        a, b = m.groups()
        c = "%s*%s" % (a, b)
    #print x, c, eval(c)
    n, m = 2, 2 # used in eval   
    return eval(c)

def convert_mol(mol, coeff):
    m = re_coeff3.match(mol)
    if m != None:
        mol, coeff = m.groups()

    coeff = convert_coeff(coeff)

    return mol, coeff

def convert_reactant_coefficients(amp):
    for r in amp.reactions:
        re = amp.reactions[r]

        new_sub_coeff = {}
        new_pro_coeff = {}

        for mol in re.substrates:
            mol, coeff = convert_mol(mol, re.sub_coeff[mol])
            new_sub_coeff[mol] = coeff

        for mol in re.products:
            mol, coeff = convert_mol(mol, re.pro_coeff[mol])
            new_pro_coeff[mol] = coeff

        #print "OLD", r, ",".join(list(re.substrates)), ",".join(list(re.products)), re.sub_coeff, re.pro_coeff
        re.substrates = set(new_sub_coeff.keys())
        re.products = set(new_pro_coeff.keys())
        re.sub_coeff = new_sub_coeff
        re.pro_coeff = new_pro_coeff
        #print "NEW", r, ",".join(list(re.substrates)), ",".join(list(re.products)), re.sub_coeff, re.pro_coeff
    #exit()

def identify_deadends(model, amp):
    mol_users = {}
    all_mols = set()
    for r in model:
        re = amp.reactions[r]
        subpro = re.substrates.union(re.products)
        all_mols.update(subpro)
        for mol in subpro:
            if mol not in mol_users:
                mol_users[mol] = set()
            mol_users[mol].add(r)

    ends = set()
    for mol in mol_users:
        if len(mol_users[mol]) <= 2: # assume both fwd and rev are in model
            ends.add(mol)
  
    return ends, all_mols

from atommap import Reaction

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

class ReactionWrapper:
    def __init__(self, reactions):
        self.reactions = reactions

def load_stoichiometry(f):
    reactions = {}
    for s in f:
        if s.startswith("#"):
            continue
        rid, balanced, origstatus, newstatus, eqn = s.strip().split("\t")
        if balanced != "True":
            #print "Skip %s" % (rid)
            continue
        lhs, rhs = eqn.split(" <=> ")
        try:
            subs, sub_coeff = parse_reactants(lhs)
            pros, pro_coeff = parse_reactants(rhs)
        except ValueError:
            continue
        r = Reaction(rid, subs, pros, sub_coeff, pro_coeff, {}, {})
        reactions[rid] = r
        rev = r.reverse()
        reactions[rev.id] = rev

    return ReactionWrapper(reactions)

def load_model(f, amp):
    model = set()
    for s in f:
        if s.startswith("#"):
            continue
        re = s.strip().split()[0].split("_")[0]  # e.g. R00001_0_rev -> R00001
        rev = "%s_rev" % (re)
        if re in amp.reactions and rev in amp.reactions:
            model.add(re)
            model.add(rev)
        else:
            #print "%s not in stoichiometry" % (re)
            pass
    return model

def load_sources(f):
    sources = {}
    for s in f:
        if s.startswith("#"):
            continue
        vals = s.strip().split()
        mol, limit = vals[2], vals[3]
        limit = float(limit)
        sources[mol] = limit
    return sources

def load_objective(f):
    objective = {}
    for s in f:
        if s.startswith("#"):
            continue
        vals = s.strip().split(None, 2)
        mol, coeff, name = vals
        coeff = float(coeff)
        if mol in objective:
            print "Error: duplicate %d in objective" % (mol)
            exit(1)
        objective[mol] = coeff
    return objective

def load_bounds(fn):
    if fn == None or fn == "-":
        return {}
    bounds = {}
    f = open(fn)
    for s in f:
        if s.startswith("#"):
            continue
        rid, lb, ub = s.strip().split()
        if lb == "-":
            lb = None
        else:
            lb = float(lb)
        if ub == "-":
            ub = None
        else:
            ub = float(ub)
        
        if lb == None:
            lb = 0.0
            revub = None
        elif lb < 0:
            revub = -lb
            lb = 0.0
        else:
            revub = 0.0
        if ub == None:
            revlb = 0.0
        elif ub < 0:
            revlb = -ub
            ub = 0.0
        else:
            revlb = 0.0

        #print rid, revlb, revub

        # Handle GLPK problem with 0,0-bounds
        if lb == 0 and ub == 0:
            ub = 1e-9
            ub = 0
        if revlb == 0 and revub == 0:
            revub = 1e-9
            revub = 0

        bounds[rid] = (lb, ub)
        bounds["%s_rev" % (rid)] = (revlb, revub)
    return bounds

def read_gibbs(f):
    f.readline()
    gibbs = {}
    for s in f:
        cid, name, dg, nh = s.strip().split(";")
        cid = int(cid)
        mol = "C%.5d" % (cid)
        gibbs[mol] = float(dg.replace(",", "."))
    return gibbs

def read_gibbs_henry(f):
    f.readline()
    gibbs = {}
    KCAL_IN_KJ = 4.19
    for s in f:
        vals = s.strip().split("\t")
        mol, dg = vals[3], vals[-2]
        if dg == "UNKNOWN":
            continue
        dg = float(dg) * KCAL_IN_KJ
        mols = mol.split("|")
        for mol in mols:
            gibbs[mol] = dg
    return gibbs

#def gibbs_to_bounds(gibbs, amp, threshold = 2225.0):
def gibbs_to_bounds(gibbs, bounds, amp, ofn, threshold = GIBBS_THRESHOLD):
    if ofn != None:
        print "Writing reaction dGs to %s..." % (ofn)
        o = open(ofn, "w")
        o.write("#Threshold=%f kJ/mol\n" % (threshold))
        o.write("#Reaction Status dG Constraint Description\n")
        #o2 = open("molecule-gibbs.txt", "w")
    missing = set()
    classes = {">" : 0, "=" : 0, "<" : 0}
    keys = amp.reactions.keys()
    keys.sort()
    ncomp = nall = 0
    nmoldg = 0
    for r in keys:
        re = amp.reactions[r]
        dg = 0.0
        dlhs = []
        localmiss = False
        
        for mol in re.sub_coeff:
            if mol in gibbs:
                dg += gibbs[mol] * re.sub_coeff[mol]
                dlhs.append("%d %s (%s)" % (re.sub_coeff[mol], mol, gibbs[mol] * re.sub_coeff[mol]))
            else:
                #print "%s not in gibbs" % (mol)
                missing.add(mol)
                localmiss = True
                dlhs.append("%d %s (?)" % (re.sub_coeff[mol], mol))
        drhs = []
        for mol in re.pro_coeff:
            if mol in gibbs:
                dg -= gibbs[mol] * re.pro_coeff[mol]
                drhs.append("%d %s (%s)" % (re.pro_coeff[mol], mol, gibbs[mol] * re.pro_coeff[mol]))
            else:
                #print "%s not in gibbs" % (mol)
                missing.add(mol)
                localmiss = True
                drhs.append("%d %s (?)" % (re.pro_coeff[mol], mol))

        if localmiss:
            rd = "="
        elif dg < -threshold:
            rd = "<"
        elif dg > threshold:
            rd = ">"
        else:
            rd = "="

        if "_rev" not in r and localmiss == False:
            classes[rd] += 1
            ncomp += 1

        if rd == "<" and r not in bounds and localmiss == False:
            bounds[r] = (0, 0)
#            print "Constraining %s" % (r)
#            pass

        if not localmiss:
            status = "complete"
        else:
            status = "missing"
        
        desc = "%s %s %s" % (" + ".join(dlhs), rd, " + ".join(drhs))

        if "_rev" not in r:
            nall += 1
            if ofn != None:
                o.write("%s\t%s\t%f\t%s\t%s\n" % (r, status, dg, rd, desc))

#    for mol in gibbs:
#        o2.write("%s\t%f\n" % (mol, gibbs[mol]))
#    for mol in missing:
#        o2.write("%s\t?\n" % (mol))

    if ofn != None:
        o.write("#%s\n" % (classes))
        o.write("#%d/%d molecules with dG (%.1f%%)\n" % (len(gibbs), len(missing) + len(gibbs), 100.0 * len(gibbs) / (len(missing) + len(gibbs))))
        o.write("#%d/%d reactions with dG (%.1f%%)\n" % (ncomp, nall, 100.0 * ncomp / nall))

#    exit()

    return bounds

def write_constraints(bounds, fn):
    f = open(fn, "w")
    keys = bounds.keys()
    keys.sort()
    f.write("#Output of %s on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))
    f.write("#Combined flux constraints from Gibbs constraints and manual fixes\n")
    for r in keys:
        if r.endswith("_rev"):
            continue
        flb, fub = bounds[r]
        rev = "%s_rev" % (r)
        if rev in bounds:
            rlb, rub = bounds[rev]
        else:
            rlb = rub = 0

        ub = fub
        lb = -rub
        #f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (r, flb, fub, rlb, rub, lb, ub))
        f.write("%s\t%s\t%s\n" % (r, lb, ub))
    f.close()

#def do_fba(kdir, adir, modelfn, sourcesfn, objectivefn, rrefn, target = None):
def do_fba(options):

    options.bounds = load_bounds(options.bounds)

    mol_names = common.read_metabolite_names(open(common.FILE_METABOLITE_NAMES))

    if options.output != None:
        o = open(options.output, "w")
        o.write("#Output of \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))
    else:
        o = sys.stdout

    err = sys.stderr
    amp = load_stoichiometry(open(options.stoichiometry))

    gibbs = None
    if options.gibbs != None:
        try:
            gibbs = read_gibbs(open(options.gibbs))
#            print "Read Gibbs energy data in Blomberg format"
        except:
            gibbs = read_gibbs_henry(open(options.gibbs))
#            print "Read Gibbs energy data in Henry format"
        options.bounds = gibbs_to_bounds(gibbs, options.bounds, amp, options.gibbs_output, options.gibbs_threshold)

    model = load_model(open(options.model), amp)

    full_model = set(amp.reactions.keys())
    if options.full:
        if options.compare:
            err.write("Warning: db compared against db (specified --full and --compare)\n")
        model = full_model
   
    rrr = set()
    for r in model:
        rrr.add(r.strip("_rev"))

    sources = load_sources(open(options.intake))

    if options.objective != None:
        objective = load_objective(open(options.objective))

    if DEADEND_SOURCES:
        deadends, all_mols = identify_deadends(model, amp)
        for mol in deadends:
            if mol not in sources and mol not in objective:
                sources[mol] = DEADEND_INTAKE

    # Modes of operation
    # 1. Test yield of individual component
    # 2. Test yield of objective
    # 3. Test yield of objective components

    if options.targets != None:
        if options.objective != None:
            err.write("Warning: using \"targets\" instead of \"objective\"\n")
        targets = options.targets.split(",")
        objective = {}
        for target in targets:
            objective[target] = -1.0

    if options.write_constraints != None:
        write_constraints(options.bounds, options.write_constraints)

    o.write("#%d reactions in model (%d bidirectional), %d sources, %d targets\n" % (len(model), len(rrr), len(sources), len(objective)))

    if options.components:
        if options.compare:
            compare_models(model, full_model, objective, sources, amp, o, mol_names, options)
        else:
            compare_models(model, None, objective, sources, amp, o, mol_names, options)
    else:
        max_yield, status = compute_yield(model, objective, sources, amp, mol_names, options)
        o.write("Yield of in model:\t%f\n" % (max_yield))
        if options.compare:
            full_max_yield, full_status = compute_yield(full_model, objective, sources, amp, mol_names, options)
            o.write("Yield in db:\t%f\n" % (full_max_yield))
            o.write("Efficiency:\t%f\n" % ((max_yield + FBA_EPS) / (full_max_yield + FBA_EPS)))
        return max_yield

    #if TEST_COMPONENTS and target == None:
    #    test_objective_components(model, objective, sources, amp)
    #elif COMPARE_FULL_VS_RECO:
    #    compare_models(full_model, model, objective, sources, amp)
    #else:
    #    max_yield, status = compute_yield(model, objective, sources, amp)
    #    print max_yield #, status

# parser = argparse.ArgumentParser(description='Perform Flux Balance Analysis on a metabolic network model.')
# #parser.add_argument('integers', metavar='N', type=int, nargs='+',
# #                   help='an integer for the accumulator')
# #parser.add_argument('--sum', dest='accumulate', action='store_const',
# #                   const=sum, default=max,
# #                   help='sum the integers (default: find the max)')
# parser.add_argument("--kegg")

# args = parser.parse_args()
# print args.kegg
# #print args.accumulate(args.integers)

VERSION = "1.0"
#DEFAULT_KEGG = "kegg"
#DEFAULT_ATOMMAP = "atom-maps-new"
#DEFAULT_STOICHIOMETRY = "balances.eqn"

USAGE_EXAMPLES = """Examples:
  %s -m model -i sources -b biomass
  %s -m model -i sources -b biomass -c
  %s -m model -i sources -b biomass -c -p
  %s -m model -i sources -t C00186
""" % (sys.argv[0], sys.argv[0], sys.argv[0], sys.argv[0])

if __name__ == "__main__":
    parser = OptionParser(version = VERSION)
    parser.add_option("-s", "--stoichiometry", 
                     help = "Stoichiometry file (default=%s)" % (common.FILE_STOICHIOMETRY),
                      metavar = "FILE")
    parser.add_option("-m", "--model",
                      help = "Metabolic network (network.reactions)", metavar = "FILE")
    parser.add_option("-i", "--intake", 
                      help = "Intake constraints", metavar = "FILE")
    parser.add_option("-t", "--targets",
                      help = "Metabolite targets")
    parser.add_option("-b", "--objective",
                      help = "Objective specification", metavar = "FILE")
    parser.add_option("-o", "--output",
                      help = "Output", metavar = "FILE")
    parser.add_option("-c", "--components", action = "store_true", default = False,
                      help = "Test objective components individually")
    parser.add_option("-f", "--full", action = "store_true", default = False,
                      help = "Use all reactions as model")
    parser.add_option("-p", "--compare", action = "store_true", default = False,
                      help = "Compare model against all reactions")
    parser.add_option("-v", "--verbose", action = "store_true", default = False, help = "Be verbose")
    parser.add_option("-l", "--fluxes", action = "store_true", default = False, help = "Print exchange flux distribution (-l -v prints all fluxes)")
    parser.add_option("-u", "--bounds", metavar = "FILE", help = "Flux lower and upper bounds")
    parser.add_option("-g", "--gibbs", metavar = "FILE", help = "Gibbs energies for metabolites")
    parser.add_option("--gibbs-threshold", help = "delta Gibbs threshold for constraining reaction directionality", type = float, default = 5.0)
    parser.add_option("--gibbs-output", help = "File to write reaction dGs to", metavar = "FILE")
    parser.add_option("--write-constraints", help = "File to write final reaction constraints to", metavar = "FILE")

    (options, args) = parser.parse_args()

    if options.stoichiometry == None or options.model == None or options.intake == None or (options.objective == None and options.targets == None):
        parser.print_help()
        print
        print USAGE_EXAMPLES
        exit()

    do_fba(options)


