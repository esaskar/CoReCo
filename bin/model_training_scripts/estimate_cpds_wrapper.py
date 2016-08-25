#!/usr/bin/env python

import sys, os, traceback

from estimate_cpds import estimate

orgf = open(sys.argv[1])   # org list
ecsdir = sys.argv[2] # where to find ec scores
modeldir = sys.argv[3]  # where to find model ECs
odir = sys.argv[4]         # output dir

for s in orgf:
    if s.startswith("#"):
        continue
    orgshort, orglong = s.strip().split()
    print orglong
    ecscoresfn = "%s/%s.ecscores" % (ecsdir, orglong)
    modelecsfn = "%s/%s" % (modeldir, orgshort)
    outfn = "%s/%s" % (odir, orgshort)
    try:
        estimate(ecscoresfn, modelecsfn, outfn)
    except Exception as e:
        print "Failed"
        traceback.print_exc(file=sys.stdout)

