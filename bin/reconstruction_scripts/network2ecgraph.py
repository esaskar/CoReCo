#!/usr/bin/env python
"""
Convert reconstruction result into EC graph.
"""

import sys, datetime
sys.path.append("../model_training_scripts/")
import common

REMOVE_COFACTORS = 1

#STOICHIOMETRY = "kegg-no-general/nooxygen.eqn"
#COFACTORS = "aux/cofactors"


STOICHIOMETRY = "../../data/Kegg/kegg-no-general/nooxygen.eqn"
COFACTORS = "../../data/Kegg/aux/cofactors"

def main(inputdir):
    ec2r, r2ec = common.read_ec_list(open(common.FILE_EC_MAP))
    reactions = common.read_stoichiometry(open(STOICHIOMETRY)).reactions
    reco = common.read_reconstruction(open("%s/%s" % (inputdir, common.NETWORK_REACTION_FILE)))
    cofactors = common.read_set(open(COFACTORS))
    o = open("%s/%s" % (inputdir, common.FILE_EC_GRAPH), "w")

    R = set()
    for r in reco:
        baser = r.split("_")[0]
        R.add(baser)

    E = {}
    for r in R:
        if r in reactions:
            re = reactions[r]
            if r not in E:
                E[r] = {}
            for m in re.substrates:
                if m not in E:
                    E[m] = {}
                E[m][r] = 1
                E[r][m] = 1
            for m in re.products:
                if m not in E:
                    E[m] = {}
                E[r][m] = 1
                E[m][r] = 1
   
    ecg = {}
    for r in E:
        if r not in r2ec:
            continue
        ec1 = r2ec[r]
        for m in E[r]:
            if m in E:
                for r2 in E[m]:
                    if r == r2:
                        continue
                    if r2 not in r2ec:
                        continue
                    ec2 = r2ec[r2]
                    for e1 in ec1:
                        if e1 not in ecg:
                            ecg[e1] = {}
                        for e2 in ec2:
                            if e2 not in ecg[e1]:
                                ecg[e1][e2] = set()
                            ecg[e1][e2].add(m)

    o.write("#Output of \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))
    cofs = list(cofactors)
    cofs.sort()
    o.write("#Cofactors: %s\n" % (",".join(cofs)))
    o.write("#EC1 EC2 SharedMetabolites AllCofactors\n")
    k1 = ecg.keys()
    k1.sort()
    for ec1 in k1:
        k2 = ecg[ec1].keys()
        k2.sort()
        for ec2 in k2:
            cofactor = 1
            for m in ecg[ec1][ec2]:
                if m not in cofactors:
                    cofactor = 0
            if cofactor and REMOVE_COFACTORS:
                continue
            o.write("%s\t%s\t%s\t%s\n" % (ec1, ec2, ",".join(ecg[ec1][ec2]), cofactor))

if __name__ == "__main__":
    main(inputdir = sys.argv[1])
