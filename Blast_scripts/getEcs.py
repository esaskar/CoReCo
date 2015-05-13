#!/usr/bin/env python

import sys, re

reEC = re.compile("DE\s+EC=([\d+\-]\.[\d+\-]\.[\d+\-]\.[\d+\-])[ ;]")

f = open(sys.argv[1])       # input: uniprot.dat
o = open(sys.argv[2], "w")  # output

idd = ac = None
ecs = []
for s in f:
    if s.startswith("ID"):
        idd = s.strip().split()[1] 
    elif s.startswith("AC"):
        ac = s.split()[1].strip().strip(";")
    elif s.startswith("DE"):
        m = reEC.match(s)
        if m:
            ecs.append(m.group(1))
    elif s.startswith("//"):
	if len(ecs) == 0:
        	e = "?"
	else:
    		e = ",".join(ecs)
        o.write("%s\t%s\t%s\n" % (ac, idd, e))
        idd = None
        ecs = []
