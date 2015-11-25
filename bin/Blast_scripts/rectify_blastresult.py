#!/usr/bin/env python

import re,sys

fn = open(sys.argv[1]) #open blast result
o = open (sys.argv[2],"w") #save output

for f in fn:
	#print f
	if f.startswith("#"):
		o.write("%s" % f)
		continue
	tmp = f.strip().split("|")

	if len(tmp) > 1:
		o.write("%s\n" % tmp[1])

	else:
		o.write("%s\n" % tmp[0])
	
