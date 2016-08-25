#!/usr/bin/env python

import re,sys,os

ln = open(sys.argv[1]) #open seq_org.list => A2R9B9  Anig
fn = open(sys.argv[2]) #open ec_file.txt => Q6GZW6  009L_FRG3G      3.6.4.-
gn = open(sys.argv[3]) #open #SeqId SeqGTGs MatchSeq MatchNID MatchGTGs MatchFrac
                       #A2P2R3  438     ?       343972  438     1.000
o = open(sys.argv[4], "w") #save result

print sys.argv[1]

orgseq={}
#save all orgseq[org]=[seq1,seq2..seqN]
for l in ln:
    if l.startswith("#"):
        continue
	seq, org = l.strip().split("\t")
	if seq not in orgseq:
		orgseq[seq] = ''
	orgseq[seq] = org
ln.close()

seqidec={}
#save all seqidec[seq]=[id,ecs]
for f in fn:
    if f.startswith("#"):
        continue
    apu = f.strip().split("\t")
    seq = apu[0]
    id = apu[1]
    if len(apu) == 3:
        ecs = apu[2]
    else:
        ecs = "?"
    #seq, id, ecs = f.strip().split("\t")
    if seq not in seqidec:
        seqidec[seq] = ''
    seqidec[seq] = ecs
fn.close()

result = {}
#extract ecs and org for each pairs of query seq VS match seq
o.write("#Org\tSeqId\tSeqGTGs\tMatchSeq\tMatchNID\tMatchGTGs\tMatchFrac\tEcs\n")
for g in gn:
	#print g.split("\t")
	if g.startswith("#"):
		continue
	try:
		qid, qgtg, mid, nid, mgtg, mfrq = g.strip().split("\t")
	except:
		print g.split("\t")
	if qid in orgseq: 
		org = orgseq[qid]
		if mid in seqidec:
			ecs = seqidec[mid]
			o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (org, qid, qgtg, mid, nid, mgtg, mfrq, ecs))
		else:
			ecs = '?'
			o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (org, qid, qgtg, mid, nid, mgtg, mfrq, ecs))
gn.close()
