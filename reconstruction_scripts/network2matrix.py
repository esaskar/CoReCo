#!/usr/bin/env python

import sys, re
sys.path.append("../model_training_scripts/")
import common


def get_reaction_id(r):
    """R00001_1_rev -> R00001"""
    return r.split("_")[0]

def write_matrix(reco, eqns, mol2name, ofn):

    omat = open("%s.matrix" % (ofn), "w")
    omol = open("%s.matrix.rows" % (ofn), "w")
    omolnames = open("%s.matrix.mols" % (ofn), "w")
    ore = open("%s.matrix.cols" % (ofn), "w")
    orenames = open("%s.matrix.reactions" % (ofn), "w")

    reactions = set()
    rid2re = {}
    for r in reco:
        rid = get_reaction_id(r)
        if not rid in eqns:
            #print "Warning: %s not in balanced equations" % (rid)
            continue
        if rid not in reactions:
            reactions.add(rid)
            rid2re[rid] = r

    mols = set()
    for rid in reactions:
        isbalanced, ostatus, nstatus, eqn = eqns[rid]
        lhs, rhs = eqn
        for mol, qty in lhs:
            mols.add(mol)
        for mol, qty in rhs:
            mols.add(mol)

    reactions = list(reactions)
    reactions.sort()
    mols = list(mols)
    mols.sort()

    N = {}
    n = len(reactions)
    for mol in mols:
        omol.write("%s\n" % (mol))
        if mol in mol2name:
            name = mol2name[mol]
        else:
            name = mol
        omolnames.write("%s\t%s\n" % (mol, name))
        N[mol] = [0 for x in range(n)]

    for ix, rid in enumerate(reactions):
        ore.write("%s\n" % (rid))

        rec = reco[rid2re[rid]]
        orenames.write("%s\t%s\n" % (rid, rec))

        isbalanced, ostatus, nstatus, eqn = eqns[rid]
        lhs, rhs = eqn
        for mol, qty in lhs:
            N[mol][ix] -= qty
        for mol, qty in rhs:
            N[mol][ix] += qty

    for mol in mols:
        row = N[mol]
        omat.write("%s\n" % ("\t".join(map(str, row))))

def convert_eqn_side(side):
    new = []
    for mol in side:
        vals = mol.split(" ")
        if len(vals) == 2:
            coeff, mol = vals
            try:
                coeff = float(coeff)
            except ValueError:
                coeff = float("NaN")
        else:
            mol = re.findall("\w\d{5}", vals[0])[0]
            #mol = vals[0]
            coeff = 1
        new.append((mol, coeff))
    return new

def convert_eqn(e):
    lhs, rhs = e.split(" <=> ")
    lhs = lhs.split(" + ")
    rhs = rhs.split(" + ")
    return (convert_eqn_side(lhs), convert_eqn_side(rhs))

def read_balanced_reactions(f):
    reactions = {}
    for s in f:
        if s.startswith("#"):
            continue
        rid, isbalanced, ostatus, nstatus, eqn = s.strip().split("\t")
        if isbalanced == "True":
            eqn = convert_eqn(eqn)
            reactions[rid] = (isbalanced, ostatus, nstatus, eqn)
    return reactions

def main(rdir, eqnfn, molfn, outfn):
    mol2name = {}
    f = open(molfn)
    for s in f:
        molid, name, name2 = s.strip().split("\t")
        mol2name[molid] = name
    print "Loading reconstruction: %s/%s" % (rdir, common.NETWORK_REACTION_FILE)
    f = open("%s/%s" % (rdir, common.NETWORK_REACTION_FILE))
    bf = open(eqnfn)
    reco = common.read_reconstruction(f)
    eqns = read_balanced_reactions(bf)
    print "%d reactions" % (len(reco))
    write_matrix(reco, eqns, mol2name, outfn)
    #sbml = convert_to_SBML(reco, eqns, mol2name)
    #libsbml.writeSBMLToFile(sbml, outfn)

if __name__ == "__main__":
    rdir = sys.argv[1]  # dir with network.reactions
    eqnfn = sys.argv[2] # balanced reactions, "balances.eqn"
    molfn = sys.argv[3]  # metabolite id -> name
    outfn = sys.argv[4] # output filename base
    main(rdir, eqnfn, molfn, outfn)
