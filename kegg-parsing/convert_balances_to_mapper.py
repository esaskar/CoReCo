#!/usr/bin/env python

import sys, datetime

f = open(sys.argv[1])  # balances.eqn
o = open(sys.argv[2], "w")  # atom mapper input

for s in f:
    if s.startswith("#"):
        continue
    rid, isbal, os, ns, eqn = s.strip().split("\t")
    if isbal == "False":
        continue
    lhs, rhs = eqn.split(" <=> ")
    lhs = map(str.split, lhs.split(" + "))
    rhs = map(str.split, rhs.split(" + "))
    subs = ""
    pros = ""
    for sub in lhs:
        if len(sub) == 2:
            coeff, mol = sub
        else:
            coeff = 1
            mol = sub[0]
        coeff = int(float(coeff))
        
        for i in range(coeff):
            subs += "%s " % (mol)
    subs = subs.rstrip()
    for pro in rhs:
        if len(pro) == 2:
            coeff, mol = pro
        else:
            coeff = 1
            mol = pro[0]
        coeff = int(float(coeff))
        for i in range(coeff):
            pros += "%s " % (mol)
    pros = pros.rstrip()
    o.write("%s %d %s p %s\n" % (rid, 1, subs, pros))
    
