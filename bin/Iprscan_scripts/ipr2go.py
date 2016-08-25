#!/usr/bin/env python

import sys, re

# open 2 files
f = open(sys.argv[1])   #input file -->.xml iprscan results
o = open(sys.argv[2],"w") # give an output file for receiving results.

iprs = {}
iprid = ''

#read input data line by line
for s in f:
	#delete all Tab and Space on the left side of each raw
	s = s.strip()
	if s.startswith("<entry ac=\"IPR"):
#<entry ac="IPR001790" desc="Ribosomal protein L10/acidic P0" name="Ribosomal_L10/acidic_P0" type="FAMILY">
	#split s by space: ['<interpro', 'id="IPR000182"', 'name="GNAT', 'domain"', 'type="Domain"', 'parent_id="IPR016181">']
		a = s.split()
	#select the second item of a and split by =
	#"IPR000182"	
		b_id = a[1].split("=")[1]
	#exlucde the ": ['', 'IPR000182', ''] and choose the second one-->get id
	#similar to get type
		iprid = b_id.split("\"")[1]
		iprs[iprid] = [iprid,[]]
		continue
        if s.find("db=\"GO\" id=\"GO:") >= 20 :
#20 is just an arbitary number
#<go-xref category="BIOLOGICAL_PROCESS" db="GO"
		a = s.split()
        #        print a[3]
	#"GO:0004635"
		b = a[3].split("=")[1]
        #        print b
	# GO:0004635
		GO = b.split("\"")[1]
        #        print  GO
		iprs[iprid][1].append(GO)
# get the keys for the dictionary data structure
keys = iprs.keys()
keys.sort()

# for ecah IPR id (here is ipr)
for ipr in keys:
	# reconbine GO term with comma 
	GO_term = ",".join(iprs[ipr][1])
	# write ipr and GO_term accordingly
	o.write("%s\t%s\n" % (ipr,GO_term))
