#!/usr/bin/env python
#cut -f 8 seqid2ec.txt | uniq | grep -E -e "[0-9]+.[0-9]+.[0-9]+.[0-9]+" | uniq | cut -d "," -f 1 | uniq | wc -l
import sys, os, re
import common
import argparse

def __main(args):
    (longtoshort, orgnames) = common.read_organisms(args.org_list)

    cutoff = 1e-10
    suffix = '.faa.IPR.final.txt'

    re_ec = re.compile("[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+")

    ddir = args.dir  # dir with org split iprscan results
    odir = args.odir # output dir

    fns = os.listdir(ddir)
    for fn in fns:
        if not fn.endswith(suffix):
                continue
        print fn    #like : name.IPR.final.txt
        f = open("%s/%s" % (ddir, fn))
        fn = fn[:-len(suffix)]   #fn like: Chaetomium_globosum
        fn = longtoshort[fn]
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

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("-d", "--dir", help = "Directory with interproscan results splitted by organism")
    p.add_argument("-o", "--odir", help = "Output directory")
    p.add_argument("--org-list", help = "org_list")
    args = p.parse_args()
    __main(args)

