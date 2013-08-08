#!/usr/bin/env python
import re,sys,os

fn = open(sys.argv[1]) #sequence from blast GTG fmt11 result: .part1
gn = open(sys.argv[2]) #starts and lens from blast GTG fmt11 result: .part2
o = open(sys.argv[3], "w") #save .result


part1 = {}
part2 = {}
for f in fn:
	if f.startswith("#"):
		continue
	Query, qname, sequence = f.strip().split("\t")
	if Query not in part1:
		part1[Query] = []
		part1[Query].append(qname)
		part1[Query].append(sequence)
fn.close()

o.write("#qname\tsname\tstarts\tlens\tqseq\n")
for g in gn:
	if g.startswith("#"):
		continue 
	Query, sname, starts, lens = g.strip().split("\t")
	try:
		o.write("%s\t%s\t%s\t%s\t%s\n" % (part1[Query][0], sname, starts, lens, part1[Query][1]))
	except:
		print "%s %s %s %s" % (Query, sname, starts, lens) 
gn.close()
