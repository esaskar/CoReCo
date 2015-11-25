#!/usr/bin/env python

import sys, os, subprocess, re
from estimate_cpds import plot

ddir = sys.argv[1]  # data dir with .pos and .neg files

fns = os.listdir(ddir)

obpos = open("%s/all.pos.blast" % (ddir), "w")
obneg = open("%s/all.neg.blast" % (ddir), "w")
ogpos = open("%s/all.pos.gtg" % (ddir), "w")
ogneg = open("%s/all.neg.gtg" % (ddir), "w")
for fn in fns:
    if fn == "all.pos.blast" or fn == "all.neg.blast" or fn == "all.pos.gtg" or fn == "all.pos.gtg":
        continue
    f = open("%s/%s" % (ddir, fn))
    if fn.endswith(".pos.blast"):
        obpos.write(f.read())
    elif fn.endswith(".neg.blast"):
        obneg.write(f.read())
    elif fn.endswith(".pos.gtg"):
        ogpos.write(f.read())
    elif fn.endswith(".neg.gtg"):
        ogneg.write(f.read())

obpos.close()
obneg.close()
ogpos.close()
ogneg.close()

plot("%s/all" % (ddir), "50 fungi")
