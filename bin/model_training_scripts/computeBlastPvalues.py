#!/usr/bin/env python

import sys, os, re, math

import common

#eps = 1e-200
approx = 1e-5

sdir = sys.argv[1]  # dir with joint blast results:[name].joint.blast
odir = sys.argv[2]  # target dir

fns = os.listdir(sdir)

for fn in fns:
    print fn
    f = open("%s/%s" % (sdir, fn))
    o = open("%s/%s.pv" % (odir, fn), "w")
    #Here I still missing One attribute 9AvgScore
    o.write("#1TargetSeq 2DBSeq 3PExt 4ECS 5BlastFwdEValue 6BlastFwdScore 7BlastRevEValue 8BlastRevScore 9AvgScore 10FwdPvalue 11RevPvalue 12JointP\n")
    c = 0
    for s in f:
        if s.startswith("#"):
            continue
	# in name.joint.blast file looks like: 8 columns
	# TargetSeq UniProtSeq PExt ECS BlastFwdEValue BlastFwdScore BlastRevEValue BlastRevScore
        seq1, seq2, pext, ecs, fwdev, fwdsc, revev, revsc = s.strip().split("\t")
        fwdp = common.blast_evalue_to_pvalue(float(fwdev))
        revp = common.blast_evalue_to_pvalue(float(revev))
        logtp = common.blast_joint_score(fwdp, revp)
        o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seq1, seq2, pext, ecs, fwdev, fwdsc, revev, revsc, fwdp, revp, logtp))

