#!/usr/bin/env python

import re, sys

fn = open(sys.argv[1]) #open blast result
o = open (sys.argv[2],"w") #save output

for s in fn:
	if s.startswith("#"):
		o.write(s)
		continue
	v = s.strip().split("\t")
	sid = v[0].split("|")
	if len(sid) > 1:
		sid = sid[1]
	else:
		sid = sid[0]
	o.write("%s\t%s\n" % (sid, "\t".join(v[1:])))
