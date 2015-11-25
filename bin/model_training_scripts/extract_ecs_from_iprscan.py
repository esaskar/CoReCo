#!/usr/bin/env python
#cut -f 8 seqid2ec.txt | uniq | grep -E -e "[0-9]+.[0-9]+.[0-9]+.[0-9]+" | uniq | cut -d "," -f 1 | uniq | wc -l
import sys, os, re

cutoff = 1e-10
suffix = '.faa.IPR.final.txt'

re_ec = re.compile("[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+")

ddir = sys.argv[1]  # dir with org split iprscan results
odir = sys.argv[2]  # output dir

fns = os.listdir(ddir)
for fn in fns:
    if not fn.endswith(".IPR.final.txt"):
            continue
    print fn    #like : name.IPR.final.txt
    f = open("%s/%s" % (ddir, fn))
    fn = fn[:-len(suffix)]   #fn like: Chaetomium_globosum
    tmp = fn[0]
    tmp1 = fn.split("_")[1][0:3]
    fn = tmp + tmp1 #short name
    print "New fn ", fn
    o = open("%s/%s" % (odir, fn), "w")
    o.write("#ECs in InterProScan results %s with cutoff = %s\n" % (fn, cutoff))
    orgecs = set()
    for s in f:
        if s.startswith("Org\t"):
            continue
        vals = s.strip().split("\t")
        evalue, ecs = vals[4], vals[7]
        ecs = ecs.split(",")
        try:
            evalue = float(evalue)
        except:
            continue
        if evalue > cutoff:
            continue
        for ec in ecs:
            if ec != "?":
                orgecs.add(ec)
    ecs = list(orgecs)
    ecs.sort()
    for ec in ecs:
        if re_ec.match(ec):
            o.write("%s\n" % (ec))
    o.close()

