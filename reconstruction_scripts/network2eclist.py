#!/usr/bin/env python
"""

"""

import sys, os, datetime
sys.path.append("../model_training_scripts/")
import common

def base(s):
    """Return reaction basename"""
    return s.split("_")[0]

def convert(reco, of):
    ec2reco = {}
    ec2r = {}
    for r in reco:
        sr = reco[r]
        ecs = sr.ec.split(",")
        for ec in ecs:
            if ec not in ec2reco:
                ec2reco[ec] = []
                ec2r[ec] = set()
            ec2reco[ec].append(sr)
            ec2r[ec].add(r)

    of.write("#Output of \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))
    of.write("#EC %s ECReactions\n" % (" ".join(common.ScoredReaction.strheader().split()[2:])))
    ecs = ec2reco.keys()
    ecs.sort()
    for ec in ecs:
        vals = "\t".join(str(ec2reco[ec][0]).split("\t")[2:])
        if len(ec2r[ec]) == 0:
            eclist = ["?_?"]
        else:
            eclist = ec2r[ec]
        of.write("%s\t%s\t%s\n" % (ec, vals, ",".join(map(base, eclist))))

def convert_dir(ddir):
    reco = common.read_reconstruction(open("%s/%s" % (ddir, common.NETWORK_REACTION_FILE)))
    convert(reco, open("%s/%s" % (ddir, common.NETWORK_ECLIST_FILE), "w"))

def convert_multiple_dirs(project_dir):
    """Convert multiple reconstruction directories under project_dir."""
    fns = os.listdir(project_dir)
    for fn in fns:
        try:
            convert_dir("%s/%s" % (project_dir, fn))
            print "Converted %s/%s" % (project_dir, fn)
        except:
            pass

def main():
    if len(sys.argv) == 3:
        ddir = sys.argv[1]  # reconstruction result directory
        ofn = sys.argv[2]   # output
        reco = common.read_reconstruction(open("%s/%s" % (ddir, common.NETWORK_REACTION_FILE)))
        convert(reco, open(ofn, "w"))
    elif len(sys.argv) == 2:
        pdir = sys.argv[1]  # project dir, reconstruction dirs as subdirs
        convert_multiple_dirs(pdir)
    else:
        print """Usage:
%s reco-dir output 
%s project-dir
""" % (sys.argv[0], sys.argv[0])

if __name__ == "__main__":
    main()
