#!/usr/bin/env python

import sys, os
from buildBlastResult import combineBlasts

RESDIR = "blast_uniprot"

f = open(sys.argv[1])  # org_list
ecfn = sys.argv[2]     # uniprot -> ec file
names = list()
for s in f:
    org, name = s.strip().split()
    names.append(name)

for name in names:
    print "Processing %s..." % (name)
    bfn1 = "%s/%s-vs-up.blast" % (RESDIR, name)
    bfn2 = "%s/up-vs-%s.blast"  % (RESDIR, name)
    bf1 = open(bfn1)
    bf2 = open(bfn2)
    ecf = open(ecfn)
    ofn = "%s/%s.joint.blast"  % (RESDIR, name)
    of = open(ofn, "w")
    combineBlasts(ecf, bf1, bf2, of)
