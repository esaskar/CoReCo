#!/usr/bin/env python
import re,sys,os

fn = open(sys.argv[1])#fmt11 blast seq
o = open(sys.argv[2], "w") # name as name.part1

Query = re.compile("local str \"(Query_.*)\"")
seq_pattern = re.compile("[A-Z]+")
seq_middle = re.compile("[A-Z]{78}")
seq_last = re.compile("[A-Z]+\"")

seq = ''
name = {}
for f in fn:
	f = f.strip()
	if Query.match(f):
		query = Query.match(f).group(1)
		if query not in name:
			name[query] = []

        #extrace query sequence name :PAS_157c00000|PAS_c159_0001
        if f.startswith("title"):
                name[query].append(f.strip().split("\"")[1].split()[0])

	#extract sequences information
	if f.startswith("seq-data ncbieaa "):
		seq = ''
		#For sequences within only one lines
		#extract sequence and save it to name[query]
		if f.endswith("\""):
			seq = f.split("\"")[1]
			name[query].append(seq)
			continue
		#This is sequence more than one line
		seq = f.split("\"")[-1]
	#middle of the sequence: SDFQEFASDQWEREADF
	if seq_middle.match(f):
		#print "find middle pattern"
		if f.endswith("\""):
			#This is the last lines
			seq = seq + f.split("\"")[0]
			name[query].append(seq)
			continue
		seq = seq + f
        #end of the sequence : OUFWODJA"
	if seq_last.match(f):
		new_seq = f.split("\"")[0]
		#print new_seq
		seq = seq + new_seq
		name[query].append(seq)
		#print name[query]

fn.close()

o.write("#Query\tqname\tqseq\n")

keys = name.keys()
for key in keys:
	#print name[key]
	try:
		o.write("%s\t%s\t%s\n" % (key, name[key][0], name[key][1]))
	except:
		print "%s %s" % (key, name[key])




