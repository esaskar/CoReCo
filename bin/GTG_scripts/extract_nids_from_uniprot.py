#!/usr/bin/env python

import re,sys,os

fn = open(sys.argv[1]) # open uniprot.fasta file
o = open(sys.argv[2], "w") # save result: nids.up

qname = set()
for f in fn:
	if f.startswith(">"):
		qname.add(f.strip().split("|")[1])
fn.close()

o.write("#nid\taccessionnumber\n")
number = 0	
for i in qname:
	number = number + 1
	o.write("%s\t%s\n" % (number, i))
o.close()					
