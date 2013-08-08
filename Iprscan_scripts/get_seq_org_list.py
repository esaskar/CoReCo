#!/usr/bin/env python
import re,sys

files = open(sys.argv[1])      #input fasta files
o = open(sys.argv[2],"w")      #output file
id = None
for f in files:
#	print "f is ",f
	if f.startswith(">"):
		id = f.split()[0].split("|")[1]
        #        print "id --->", id
		o.write("%s\tPsti\n" % (id))

		
