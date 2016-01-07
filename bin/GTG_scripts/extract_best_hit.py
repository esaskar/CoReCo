#!/usr/bin/env python

import re,sys,os
#Ustilago_maydis.faa-final-input.txt
fn = open(sys.argv[1]) # read -final-input.txt
o = open(sys.argv[2], "w")  # save result

besthit = {}
#divide file into group divided by qname
for f in fn:
	if f.startswith("#"):
		continue
	#fn, sn, q.name, q.start, q.end, rep-level, nid, s.start, s.end, bitscore, length, %identity, detial pair information(start_1, start_2, lens), q.seq
	fns, sn, qname, qstart, qend, rep, nid, sstart, send, bitscore, length, identity, start1, start2, lens, qseq = f.strip().split("\t")
	#print "%s %s" % (qname, nid)
	if qname not in besthit:
		i = 0
		besthit[qname] = {}
	if len(besthit[qname]) == 0:
		besthit[qname][i] = [fns, sn, qname, qstart, qend, rep, nid, sstart, send, bitscore, length, identity, start1, start2, lens, qseq]
		i = i + 1
		continue
	if qname in besthit:
		tmp = [fns, sn, qname, qstart, qend, rep, nid, sstart, send, bitscore, length, identity, start1, start2, lens, qseq]
		besthit[qname][i] = tmp
		i = i + 1
		#print qname
fn.close()

qnames = besthit.keys()
print len(qnames)
qnames.sort()

maxscore = {}
#choose the maxscore of each query sequence
for q in qnames:
	if q not in maxscore:
		maxscore[q] = 0
	for k in besthit[q].keys():
		if float(besthit[q][k][9]) > float(maxscore[q]):
			maxscore[q] = besthit[q][k][9] 

o.write("#fn, sn, q.name, q.start, q.end, rep-level, nid, s.start, s.end, bitscore, length, %identity, detial pair information(start_1, start_2, lens), q.seq\n")
for key in qnames:
	for k in besthit[key].keys():
		if float(besthit[key][k][9]) == float(maxscore[key]):
			tmp = "\t".join(besthit[key][k])
			o.write("%s\n" % tmp)
o.close()










	

