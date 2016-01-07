#!/usr/bin/env python

import os, re, sys

gn = open(sys.argv[1]) # blast output format 6 with qseq
fn = open(sys.argv[2]) # blast output (from 11 to my form): containing qname sname start lens
o = open(sys.argv[3], "w") # save results

result_fmt6 = {}
count_fmt6 = 0
#./blast_formatter -archive test11 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq'
# saving attritubtes of fn1 into lists
# This is for format 6
for g in gn:
        #qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = f.strip().split("\t")
        if g.startswith("#"):
                continue

        count_fmt6 = count_fmt6 +1
        result_fmt6[count_fmt6] = g.strip()
#print  "This is count_fmt6 %s " % count_fmt6

#define a dictionary to save result from fn file
#save information from fn in result; result[qname][sname] = [starts, lens, seq]
result = {}

for f in fn:
	if f.startswith("#"):
		continue
	qname, sname, starts, lens, qseq = f.strip().split("\t")
	if qname not in result:
		result[qname] = {}
	if sname not in result[qname]:
		result[qname][sname] = {}
	#same query and subject seq can be paired at different position, so I add this.
	if starts not in result[qname][sname]:
		result[qname][sname][starts] = [lens, qseq]
fn.close()	


#still missing nid becasue no blast result
o.write("#fn, sn, q.name, q.start, q.end, rep-level, nid, s.start, s.end, bitscore, length, %identity, detial pair information(start_1, start_2, lens), q.seq\n")

nid_pattern = re.compile("[0-9]+")
starts_1 = ''
starts_2 = ''
a = []
t = 0
for key in range(1, count_fmt6+1):
        #print "key is %s" % key
        #print "result_fmt6[key] is %s" % result_fmt6[key]
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = result_fmt6[key].split("\t")
        #nid is the group of number before the first | of sseqid        
        nid = sseqid.strip().split("|")[0]
        if not nid_pattern.match(str(nid)):
                nid = "?"
        #print "starts[key] is %s" % starts[key]  
	if qseqid in result.keys() and sseqid in result[qseqid].keys():
		keys = result[qseqid][sseqid].keys()
		
		for s in keys:     
        		tmp = s.strip().split()
        		for n in range(len(tmp)):
                		if len(starts_1) == 0:
					if n%2 == 0:
                	                	starts_1 = tmp[n]
                		if len(starts_2) == 0:
                		        if n%2 != 0:
						starts_2 = tmp[n]
                		else:
                		        if n%2 == 0:
                	        	        starts_1 = starts_1 + ' ' + tmp[n]
                	        	if n%2 != 0:
                	                	starts_2 = starts_2 + ' ' + tmp[n]
			#print "starts_1 is %s\nstarts_2 is %s\nlens is %s\n" % (starts_1, starts_2, result[qseqid][sseqid][1])
			starts_1 = starts_1.split()
			starts_2 = starts_2.split()
			lens = result[qseqid][sseqid][s][0].split()
			#print "a is %s\nstarts_1 is %s\nstarts_2 is %s\nlens is %s\n" % (a, starts_1, starts_2, lens)
			#Delete all itmes of -1 from starts_1 and starts_2 and the same index of lens
			#record index of -1 in starts_1
			#following if is used to distinguish differnt matches of differnt starts for the same query and subject seq
			#print "qname ", qseqid, "sname ", sseqid, "qstart", qstart, "qend ", qend, "starts_1[0] ", (int(starts_1[0])+1), "starts_1[-1]+lens[-1] ", int(starts_1[-1])+int(lens[-1])
			if (int(starts_1[0])+1) == int(qstart) and (int (starts_1[-1])+int(lens[-1])) == int(qend) and (int(starts_2[0])+1) == int(sstart) and (int(starts_2[-1])+int(lens[-1])) == int(send):
				for i in range(len(starts_1)):
					if starts_1[i] == '-1':
						a.append(i)
				#delet all items in starts_1 and the same index of starts_2 and lens if a > 1
				#print "a is %s\nstarts_1 is %s\nstarts_2 is %s\nlens is %s\n" % (a, starts_1, starts_2, lens)
				if len(a) > 0:
					for i in range(len(a)):
						try:
							starts_1[a[i]] = []
							starts_2[a[i]] = []
							lens[a[i]] = []
						except:
							print "i %d\na[i] %d\nstarts_1 %s\nstarts_2 %s\nlens %s\nstarts_1[a[i]] %s\nstarts_2[a[i]] %s" % (i, a[i], starts_1, starts_2, lens, starts_1[a[i]], starts_2[a[i]])
					a = []
				#delete all items in starts_2 and the same index of starts_1 and lens if a > 1
				for j in range(len(starts_2)):
					if starts_2[j] == '-1':
						a.append(j)
				#print "a is %s\nstarts_1 is %s\nstarts_2 is %s\nlens is %s\n" % (a, starts_1, starts_2, lens)
				if len(a) > 0:
					for j in range(len(a)):
					#	print a[j]
						try:
							starts_1[a[j]] = []
							starts_2[a[j]] = []
							lens[a[j]] = []
						except:
							print "j %d\na[j] %d\nstarts_1 %s\nstarts_2 %s\nlens %s\nstarts_1[a[j]] %s\nstarts_2[a[j]] %s" % (j, a[j], starts_1, starts_2, lens, starts_1[a[j]], starts_2[a[j]])

					a = []
				starts_1 = filter(None, starts_1)
				starts_2 = filter(None, starts_2)
				lens = filter(None, lens)
			#	print starts_1
			#	print starts_2
			#	print lens
				starts_1 = " ".join(starts_1)
				starts_2 = " ".join(starts_2)
				lens = " ".join(lens)
        			o.write("1\t%s\t%s\t%s\t%s\tnrdb40\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (key, qseqid, qstart, qend, nid, sstart, send, bitscore, length, pident, (starts_1+'\t'+starts_2), lens + '\t' + result[qseqid][sseqid][s][1]))        
			starts_1 = ''
			starts_2 = ''

