#!/usr/bin/python
#
# makes a python dictionary out of the "reaction" file in kegg ligand
#

import re



REACTIONFILE = "reaction"
R = {}


f = open(REACTIONFILE)
lines = map(str.strip, f.readlines())
f.close()


r = None


for i,line in enumerate(lines):
	prefix = line[0:13].strip()
	data = line[13:].strip()
	
	if prefix == "ENTRY":
		r = data.split()[0]
		R[r] = {}
	elif prefix:
		R[r][prefix]











