#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
"""
Balance KEGG reactions when possible by finding reaction coefficients such that left and right hand sides of reaction equation have equal number of the same atom type and the sum of coefficients is minimized.

We allow certain molecules to be freely added to either side in reaction, see AUX_MOLS.

Requires GLPK and PyGLPK.

Author: Esa Pitkänen (www.cs.helsinki.fi/u/epitkane)
"""

import sys, re, optparse
import glpk

#NO_BALANCING = 0               # 1 = just write unmodified stoichiometries to output
#IGNORE_MISSING_FORMULAE = 0    # 1 = balance reactions with missing formulae; 0 = do not balance reaction if any formula is missing
#RETAIN_MISSING = 0             # 1 = keep reactions with missing formulae unmodified

# Maximum reactant coefficient allowed
MAX_COEFF = 15

# Molecules that are freely available in balancing
AUX_MOLS = {"C00001" : MAX_COEFF,  # H2O
            "C00080" : MAX_COEFF,  # H+
            "C01353" : 1          # Carbonic acid (C1)
#            "C00007" : MAX_COEFF,  # Oxygen
#            "A00001" : MAX_COEFF   # Wildcard *
}  

# e.g., "C", "C3", "Cu", "*"
re_elem = re.compile("([A-Z\*][a-z]*)(\d*)")

def parse_formula(s):
    """Parse a compound formula in KEGG. Does not handle generic compounds such as C10H17O8PR2(C5H8O5PR)n very well."""
    # assume n=1 in formula
    # C4(C7)n -> C4C7 = C11
    #print "***", s
    s = s.translate(None, ". ()n")
    elems = {}
    for elem, qty in re_elem.findall(s):
        if qty == "":
            qty = 1
        else:
            qty = int(qty)
        if elem not in elems:
            elems[elem] = 0
        elems[elem] += qty
    return elems

def read_compound(f):
    """Parse KEGG compound file."""
    formulae = {}
    entry = None
    for s in f:
        if s.startswith("ENTRY"):
            entry = s.strip().split()[1]
            formula = None
        elif s.startswith("FORMULA"):
            formula = s.strip().split()[1]
        elif s.startswith("///"):
            if formula != None:
                formulae[entry] = parse_formula(formula)
    return formulae

re_molparen = re.compile("([\w\d]+)(\(.+\))?")
re_mn = re.compile("(\d+)([nm])")

def side_balance(reactants, bal, factor, formulae, options):
    #print "side_balance", reactants
    for mol in reactants:

        vals = mol.split()
        if len(vals) == 1:
            qty, mol = 1, vals[0]
        else:
            qty, mol = vals

        vals = re_molparen.findall(mol)[0]
        mol, qty2 = vals
        if qty2 != "":
            qty = qty2
        qty = str(qty)
        m = re_mn.match(qty)
        if m != None:
            a, b = m.groups()
            qty = "%s*%s" % (a, b)
        n, m = 1, 1
        qty = eval(str(qty))

        if mol not in formulae:
            if not options.ignore_missing:
                return False
            continue

        for elem in formulae[mol]:
            if elem not in bal:
                bal[elem] = 0
            bal[elem] += factor * qty * formulae[mol][elem]

    return True

def convert_nonnumeric_equation(eqn):
    return eqn.replace("n", "1").replace("m", "1") # :p

def convert_quantity(qty):
    try:
        qty = int(qty)
    except:
        try:
            qty = int(eval(qty))
        except:
            qty = 1
    return qty

def reactants_to_atoms(reactants, mult, formulae, aux_mols, options):
    stuff = {}
    molqty = {}
    missing = False
    for mol in reactants:
        vals = mol.split()
        if len(vals) == 1:
            qty, mol = 1, vals[0]
        else:
            qty, mol = vals
        if mol not in molqty:
            molqty[mol] = 0
        qty = convert_quantity(qty)
        molqty[mol] += int(qty)
        if mol not in formulae:
            missing = True
            if options.ignore_missing:
                continue
            else:
                return False, {}, [], True
        stuff[mol] = {}
        for elem in formulae[mol]:
            stuff[mol][elem] = mult * formulae[mol][elem]
    for mol in aux_mols:
        if mol not in stuff:
            stuff[mol] = {}
            for elem in formulae[mol]:
                stuff[mol][elem] = mult * formulae[mol][elem]

    return True, stuff, molqty, missing

def glpkout(s):
    print "glpk: %s" % (s)


def balance_reaction(eqn, formulae, options):
    lhs, rhs = eqn.split(" <=> ")
    lhs = lhs.split(" + ")
    rhs = rhs.split(" + ")
    bal = {}

    # Find out the atom composition of reactants
    ok, lhs_atoms, lhs_qty, lhs_missing = reactants_to_atoms(lhs, -1, formulae, AUX_MOLS, options)
    if not ok:
        return None, "missing"  # cannot balance: missing reactant formula
    ok, rhs_atoms, rhs_qty, rhs_missing = reactants_to_atoms(rhs, 1, formulae, AUX_MOLS, options)
    if not ok:
        return None, "missing"

    if options.retain_missing and (lhs_missing or rhs_missing):
        return eqn, "missing"

    reactant2ix = {}
    ix2reactant = {}
    atom2ix = {}
    ix = 0
    aix = 0
    N = {}

    # Prepare a balancing matrix N
    for mol in lhs_atoms:
        reactant2ix["L-%s" % (mol)] = ix
        ix2reactant[ix] = "L-%s" % (mol)
        for atom in lhs_atoms[mol]:
            if atom not in atom2ix:
                atom2ix[atom] = aix
                aix += 1
            N[(ix, atom2ix[atom])] = lhs_atoms[mol][atom]
        ix += 1
    for mol in rhs_atoms:
        reactant2ix["R-%s" % (mol)] = ix
        ix2reactant[ix] = "R-%s" % (mol)
        for atom in rhs_atoms[mol]:
            if atom not in atom2ix:
                atom2ix[atom] = aix
                aix += 1
            N[(ix, atom2ix[atom])] = rhs_atoms[mol][atom]
        ix += 1

    lp = glpk.LPX()
    lp.name = "b"
    lp.obj.maximize = False
    lp.rows.add(aix)
    lp.cols.add(ix)
    for i in range(aix):
        lp.rows[i].bounds = 0, 0  # atom balances
    for i in range(ix):
        side, mol = ix2reactant[i].split("-")
        # if mol in AUX_MOLS:
        #     lp.cols[i].bounds = 0, AUX_MOLS[mol]
        # else:
        #     if side == "L":
        #         qty = lhs_qty[mol]
        #     else:
        #         qty = rhs_qty[mol]
        #     lp.cols[i].bounds = qty, max(qty, MAX_COEFF) # "1, None" crashes glpk intopt
        if mol in lhs_qty and side == "L":
            qty = lhs_qty[mol]
        elif mol in rhs_qty and side == "R":
            qty = rhs_qty[mol]
        else:
            qty = 0
        if mol in AUX_MOLS:
            limit = AUX_MOLS[mol]
        else:
            limit = MAX_COEFF
        lp.cols[i].bounds = qty, max(qty, limit) # "1, None" crashes glpk intopt
        lp.cols[i].name = ix2reactant[i]

    lp.obj[:] = [1 for x in range(ix)]
    M = []
    for u, v in N:
        M.append((v, u, N[(u, v)]))
    lp.matrix = M
    
    glpk.env.term_on = False
    #glpk.env.term_hook = glpkout

    # glpk integer solver requires us to solve linear problem first
    retval = lp.simplex()
    assert retval == None 
    if lp.status != 'opt': 
        return None, "no-solution"

    for col in lp.cols:
        col.kind = int

    lp.intopt()

    # Convert solution to equation format
    blhs = []
    brhs = []
    all_zeros = True
    for c in lp.cols:
        if c.primal != 0:
            all_zeros = False
            if c.name[0] == "L":
                blhs.append("%s %s" % (c.primal, c.name.split("-")[1]))
            else:
                brhs.append("%s %s" % (c.primal, c.name.split("-")[1]))

    if all_zeros:
        return None, "no-solution"

    beq = " + ".join(blhs)
    for mol in lhs_qty:
        if mol not in lhs_atoms:
            beq += " + %d %s" % (lhs_qty[mol], mol)
    beq += " <=> " + " + ".join(brhs)
    for mol in rhs_qty:
        if mol not in rhs_atoms:
            beq += " + %d %s" % (rhs_qty[mol], mol)
    return beq, "balanced"

def get_balance_code(balances):
    n = len(balances.keys())
    if n == 0:
        st = "missing"
    elif "H" in balances and sum(balances.values()) == balances["H"]:
        st = "proton"
    else:
        st = "imbalance"
    return st

def check_balance(eqn, formulae, options):
    lhs, rhs = eqn.split(" <=> ")
    lhs = lhs.split(" + ")
    rhs = rhs.split(" + ")
    bal = {}

    if not side_balance(lhs, bal, -1, formulae, options):
        return False, {}, "missing"
    if not side_balance(rhs, bal, 1, formulae, options):
        return False, {}, "missing"

    balanced = True
    for mol in bal:
        if bal[mol] != 0:
            balanced = False
            break

    if not balanced:
        st = get_balance_code(bal)
    else:
        st = "balanced"

    return balanced, bal, st

def check_reactions(f, formulae, options):
    eqnon = 0

    of = open("%s.status" % (options.output), "w")
    of2 = open("%s.eqn" % (options.output), "w")
    import datetime
    msg = "#Output of \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now())

    msg += "#no-balancing = %s\n#ignore-missing-formulae = %s\n#retain-missing = %s\n" % (options.no_balancing, options.ignore_missing, options.retain_missing)

    of.write(msg)
    of2.write(msg)
    of2.write("#Reaction Balanced OrigStatus NewStatus Equation\n")

    status = {"balanced" : 0,
              "missing" : 0,
              "imbalance" : 0,
              "proton" : 0,
              "total" : 0}

    new_status = {"balanced" : 0,
              "missing" : 0,
              "imbalance" : 0,
              "proton" : 0,
              "total" : 0,
              "cannot-balance" : 0,
              "already-balanced" : 0}

    c = 0
    for s in f:
        if s.startswith("ENTRY"):
            entry = s.strip().split()[1]
        elif s.startswith("EQUATION"):
            eqn = s[12:].strip()
            eqnon = 1
        elif s.startswith("///"):
            c += 1

            if options.no_balancing:
                of.write("   %s (%s)\n" % (eqn, "unmodified"))
                of2.write("%s\t%s\t%s\t%s\t%s\n" % (entry, True, "unmodified", "balanced", eqn))
                continue

            success = False
            eqn = convert_nonnumeric_equation(eqn) # n,m -> 1
            orig_balanced, orig_balances, orig_code = check_balance(eqn, formulae, options)
            # if entry == "R00002":
            #     print eqn
            #     print orig_balanced
            #     print orig_balances
            #     print orig_code

            status["total"] += 1

            new_eqn = new_balanced = None
            new_code = "already-balanced"

            if not orig_balanced:
                new_eqn, new_balanced = balance_reaction(eqn, formulae, options)
                if new_eqn == None:
                    new_code = "cannot-balance"
                # if entry == "R00002":
                #     print "============="
                #     print new_eqn
                #     print new_balanced
            else:
                success = True

            of.write("%s\t%s\n" % (entry, orig_code))
            of.write("   %s\n" % (eqn))
            of.write("    balances: %s\n" % (orig_balances))
            of.write("    balancing: %s\n" % (new_balanced))
            if new_eqn != None:
                new_balanced, new_balances, new_code = check_balance(new_eqn, formulae, options)
                if new_balanced:
                    success = True
                # if entry == "R00002":
                #     print new_balanced
                #     print new_balances
                #     print new_code
                #     exit()

                of.write("   %s (%s)\n" % (new_eqn, new_code))
                of2.write("%s\t%s\t%s\t%s\t%s\n" % (entry, success, orig_code, new_code, new_eqn))
            else:
                of2.write("%s\t%s\t%s\t%s\t%s\n" % (entry, success, orig_code, new_code, eqn))

            status[orig_code] += 1
            new_status[new_code] += 1

        elif s[0] == " " and eqnon:
            eqn += " " + s.strip()
        else:
            eqnon = 0
   
    of.write("#Originals: %s\n" % (status))
    of.write("#Balanced: %s\n" % (new_status))

def read_formulae(f):
    formulae = {}
    for s in f:
        if s.startswith("#"):
            continue
        vals = s.strip().split("\t")
        if len(vals) == 3:
            mol, ss, form = vals
        elif len(vals) == 2:
            mol, ss = vals
        else:
            assert(0)
        if ss == "?":
            continue
        ss = ss.split(",")
        formulae[mol] = {}
        for val in ss:
            atom, qty = val.split(":")
            qty = int(qty)
            formulae[mol][atom] = qty
    return formulae

def add_dummy_formulae(formulae):
    #formulae["ANYMOL"] = {"*" : 1}
    pass
        
def main(options):
    #formulae = read_compound(open("%s/compound" % (kdir)))
    #formulae = read_formulae(open("aux/kegg-metabolite-formulae2"))

    formulae = read_formulae(open("../../data/Kegg/aux/kegg-metabolite-formulae"))
    #formulae = read_formulae(open("%s/aux/kegg-metabolite-formulae"%(options.kegg))

    add_dummy_formulae(formulae)
    check_reactions(open("%s/reaction" % (options.kegg)), formulae, options)

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-k", "--kegg", 
                      help = "KEGG directory", metavar = "DIR")
    parser.add_option("-o", "--output",
                      help = "Output file", metavar = "FILE")
    parser.add_option("--no-balancing", action = "store_true", default = False,
                      help = "Write unmodified reactions to output")
    parser.add_option("--ignore-missing", action = "store_true", default = False,
                      help = "Balance reactions with missing formulae")
    parser.add_option("--retain-missing", action = "store_true", default = False,
                      help = "Keep reactions with missing formulae unmodified")
    
    (options, args) = parser.parse_args()
    if not options.kegg or not options.output:
        parser.print_help()
        exit(2)
    main(options)
