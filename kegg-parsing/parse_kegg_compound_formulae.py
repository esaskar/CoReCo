#!/usr/bin/env python

import sys, os, datetime

ddir = sys.argv[1]  # mol dir
ofn = sys.argv[2]   # output file

o = open(ofn, "w")
o.write("#Output of \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))

fns = os.listdir(ddir)
for fn in fns:
    atoms = {}
    f = open("%s/%s" % (ddir, fn))
    f.readline()
    f.readline()
    f.readline()
    header = f.readline()
    na = int(header[0:3])
    for i in range(na):
        vals = f.readline().strip().split()
        atom = vals[3]
        if atom not in atoms:
            atoms[atom] = 1
        else:
            atoms[atom] += 1
    k = atoms.keys()
    k.sort()
    ss = ""
    for atom in k:
        ss += "%s:%d," % (atom, atoms[atom])
    o.write("%s\t%s\n" % (fn[:-4], ss.rstrip(",")))
