#!/usr/bin/env python

import sys

kdir = sys.argv[1]
o = open(sys.argv[2], "w")   # EC -> pathway
o2 = open(sys.argv[3], "w")  # pathway -> name
o3 = open(sys.argv[4], "w")  # reaction -> name
f = open("%s/reaction" % (kdir))
pathon = 0
ec2path = {}
r2path = {}
path2name = {}
for s in f:
    if s[0] != " ":
        pathon = 0
    if s.startswith("ENTRY"):
        entry = s.strip().split()[1]
        path = set()
        ec = None
    elif s.startswith("PATHWAY"):
        pathon = 1
        vals = s[12:].strip().split(None, 1)
        if len(vals) < 2:
            continue
        pid, pname = vals
        path.add(pid)
        if pid not in path2name:
            path2name[pid] = pname
    elif s[0] == " " and pathon:
        vals = s[12:].strip().split(None, 1)
        if len(vals) < 2:
            continue
        pid, pname = vals
        path.add(pid)
        if pid not in path2name:
            path2name[pid] = pname
    elif s.startswith("ENZYME"):
        ec = s.strip().split()[1]
        if ec not in ec2path:
            ec2path[ec] = set()
    elif s.startswith("///"):
        r2path[entry] = set()
        for pid in path:
            r2path[entry].add(pid)
            if ec:
                ec2path[ec].add(pid)

keys = ec2path.keys()
keys.sort()
for ec in keys:
    o.write("%s\t%s\n" % (ec, ",".join(list(ec2path[ec]))))

keys = path2name.keys()
keys.sort()
for pid in keys:
    o2.write("%s\t%s\n" % (pid, path2name[pid]))

keys = r2path.keys()
keys.sort()
for r in keys:
    o3.write("%s\t%s\n" % (r, ",".join(list(r2path[r]))))
