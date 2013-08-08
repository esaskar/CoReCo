#!/usr/bin/env python

import sys, os

def main(ddir, ofn):
    fns = os.listdir(ddir)
    o = open(ofn, "w")
    fns.sort()
    data = {}
    mkeys = set()
    for fn in fns:
        f = open("%s/%s" % (ddir, fn))
        data[fn] = {}
        for s in f:
            if s.startswith("#"):
                continue
            vals = s.strip().split("\t")
            mol, mname, val = vals[0:3]
            data[fn][mname] = val
            mkeys.add(mname)

    mkeys = list(mkeys)
    mkeys.sort()
    o.write("#")
    for mk in mkeys:
        o.write("%s " % (mk.replace(" ", "_")))
    o.write("\n")
    for fn in fns:
        if "." not in fn:
            base = fn
        else:
            base = ".".join(fn.split(".")[0:-1])
        o.write("%s" % (base))
        for mk in mkeys:
            if mk not in data[fn]:
                o.write("\t?")
            else:
                o.write("\t%s" % (data[fn][mk]))
        o.write("\n")

if __name__ == "__main__":
    ddir = sys.argv[1] # directory with fba.py result files <species>.txt
    ofn = sys.argv[2]  # output
    main(ddir, ofn)
