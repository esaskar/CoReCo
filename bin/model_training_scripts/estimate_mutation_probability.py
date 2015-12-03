#!/usr/bin/env python

import sys, os

import common, fitch, tree

def root_tree(tree, root):
    children = {}
    Q = [root]
    while len(Q) > 0:
        q = Q.pop(0)
        children[q] = {}
        for u in tree[q]:
            if u in T:
                continue
            children[q][u] = 1
            parent
            Q.append(u)
    return parents, children

def estimate(mdir, species_tree_fn, odir):
    of = open("%s/%s" % (odir, common.TREE_CPD_FILE), "w")
    models, all_ecs = common.read_models(mdir)
    t = common.read_tree(species_tree_fn)   
    parents, children = tree.findParentsAndChildren(t, "N1")
    snadd, sndel, snnoadd, snnodel = 0, 0, 0, 0
    all_mutations = {}
    for i in range(500):
        nadd, ndel, nnoadd, nnodel, mutations = fitch.parsimony(parents, children, models, all_ecs)
        snadd += nadd
        sndel += ndel
        snnoadd += nnoadd
        snnodel += nnodel
        for key in mutations:
            if key not in all_mutations:
                all_mutations[key] = [0, 0, 0, 0]
            all_mutations[key][0] += mutations[key][0] + 1
            all_mutations[key][1] += mutations[key][1] + 1
            all_mutations[key][2] += mutations[key][2] + 1
            all_mutations[key][3] += mutations[key][3] + 1

    #print all_mutations

    keys = all_mutations.keys()
    keys.sort()
    for key in keys:
        nadd, ndel, nnadd, nndel = all_mutations[key]
        padd = 1.0 * nadd / (nadd + nnadd)
        pdel = 1.0 * ndel / (ndel + nndel)
        u, v = key.split("-")
        of.write("%s\t%s\t%s\t%s\n" % (u, v, padd, pdel))
        #print key, padd, pdel, nadd, ndel, nnadd, nndel

    #print "p-add\t%s\t%s" % (1.0 * snadd / (snadd + snnoadd), snadd)
    #print "p-del\t%s\t%s" % (1.0 * sndel / (sndel + snnodel), sndel)
    #print "p-no-add\t%s\t%s" % (1.0 * snnoadd / (snadd + snnoadd), snnoadd)
    #print "p-no-del\t%s\t%s" % (1.0 * snnodel / (sndel + snnodel), snnodel)
    of.close()

if __name__ == "__main__":
    modeldir = sys.argv[1]  # dir with models
    species_tree_fn = sys.argv[2]      # species tree
    odir = sys.argv[3] # output dir
    estimate(modeldir, species_tree_fn, odir)
