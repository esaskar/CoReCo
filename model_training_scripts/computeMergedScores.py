#!/usr/bin/env python
#
# Merge BLAST and GTG scores into .rawscores
#

import sys, os, re

f = open(sys.argv[1]) # organism list, e.g., org_list
blastdir = sys.argv[2] # the dir contains name.joint.blast.pv files of 50 fungis
gtgdir = sys.argv[3] # gtg score dir I guess... right now I don't have it, change the number ahead
outdir = sys.argv[4] # save the result

BASELINE_JOINTP = 1 
BASELINE_GTGSCORE = 0.1
#user can define a cutoff for blast and gtg by themself
if len(sys.argv) > 6:
    blastscorecutoff = float(sys.argv[5])
    gtgscorecutoff = float(sys.argv[6])
else:
    blastscorecutoff = None
    gtgscorecutoff = None
longtoshort = {}
#in org_list: data looks like: Pichia  pichia_stipistic
#examples:
#Anig Aspergillus_niger
#Tree Trichoderma_reesei
#Pgra Puccinia_graminis
for s in f:
    if s.startswith("#"):
        continue
    shortn, longn = s.strip().split(" ")
    longtoshort[longn] = shortn

if len(longtoshort) == 0:
    sys.stderr.write("Warning: no organisms in %s\n" % (sys.argv[1]))

# save all of the -join.blast file as a list from blastdir in fns
fns = os.listdir(blastdir)
fns.sort()
#Here the name.joint.blast.pv file is the result of name.joint.blast from 
#computeBlastPvalues.py script.
# "#"the original file
reblastfn = re.compile("(.+)\.joint\.blast\.pv")


i = 0
for fn in fns:
    i += 1
    orgname = reblastfn.findall(fn) # orgname save the long name of each species
    orgname = orgname[0]
    print "orgname--->", orgname
    
    try:
        short = longtoshort[orgname]
	
    except Exception:
	
	short = None

    if (not short is None) and ( not os.path.exists("%s/%s.rawscores" % (outdir, orgname)) ):	
	#show following sentence on terminal:1/total number: Merging .... of orgname
	 sys.stdout.write("%d/%d: Merging scores of %s: " % (i, len(fns), orgname))
	 sys.stdout.flush()
	#f: saving the file: name.joint.blast.pv
	 f = open("%s/%s" % (blastdir, fn))
	#give the output [name].rawscores; bulid if not exist
	#give the output [longname].kde; bulid if not exist
	 o = open("%s/%s.rawscores" % (outdir, orgname), "w")  # joint pvalues and gtg scores combined
	#I have no idea what is kde score
	 o2 = open("%s/%s.kde" % (outdir, orgname), "w") # kde scores
	 o2.write("#Seq EC MaxPValue MaxGTGScore KernelDensity\n")
	 #save [seq1][seq2]=jointp
	 blastpvalues = {}
	 #dbseqtoec: find all ecs for seq2
	 dbseqtoec = {}
	 #find all seq2 for seq1
	 allseqs = {} # seq -> seq -> 1
     
	 # read blast results & pvalues for a single species
	 f.readline() # header
	 sys.stdout.write("BLAST ")
	 sys.stdout.flush()
	 print "* Reading BLAST..."
	 for s in f:
	     #print "s-->",s
	     if s.startswith("#"):
		 continue
	     #seq1, seq2, pext, ecs, evf, scf, evr, scr, avgsc, pvf, pvr, jointp = s.strip().split("\t")
	     #acturally I only contain the first 8 attributes.So I "#" the original one
	     seq1, seq2, pext, ecs, evf, scf, evr, scr, pvf, pvr, jointp = s.strip().split("\t")
	     #seq1, seq2, pext, ecs, evf, scf, evr, scr 
	     #seq1, seq2, pext, ecs, evf, scf, evr, scr = s.strip().split("\t")
	     if seq2 not in dbseqtoec:
		 #findout what is set():to be unique
		 dbseqtoec[seq2] = set()
	     dbseqtoec[seq2].add(ecs)
	     #allseqs[seq1]=[seq2_1,seq2_2.....] for example
	     if seq1 not in blastpvalues:
		 blastpvalues[seq1] = {}
		 allseqs[seq1] = set()
	     #allseqs[seq1] looks like: set([seq2]) after .add()
	     allseqs[seq1].add(seq2)
	     #select the joint.blast result over the minimum blastscorecutoff.
	     #blastpvalues[seq1][seq2]= jointp for example
	     jointp = float(jointp)
	     if blastscorecutoff == None or jointp < blastscorecutoff:
		 blastpvalues[seq1][seq2] = jointp
	     else:
		 blastpvalues[seq1][seq2] = 0
     
	 # read gtg results for the same species
	 # and try to find those miss by blast result but recognise y GTG feature
	 # I need to change the name into short.gtg.knn
	 gtgf = open("%s/%s.gtg.knn" % (gtgdir, short))
	 gtgf.readline() # header
	 gtgscores = {}
	 print "* Reading GTG..."
	 sys.stdout.write("GTG ")
	 sys.stdout.flush()
	
	 for s in gtgf:
	     if s.startswith("#"):
		 continue
	     org, seq1, seqgtg, seq2, matchnid, matchgtg, matchfrac, ecs = s.strip().split("\t")
	     if seq1 not in allseqs:
		 allseqs[seq1] = set()
	     allseqs[seq1].add(seq2)
	     if seq1 not in gtgscores:
		 gtgscores[seq1] = {}
	     if seq2 not in dbseqtoec:
		 dbseqtoec[seq2] = set()
	     dbseqtoec[seq2].add(ecs)
     
	     matchfrac = float(matchfrac)
	     if gtgscorecutoff == None or matchfrac < gtgscorecutoff:
		 gtgscores[seq1][seq2] = matchfrac
	     else:
		 gtgscores[seq1][seq2] = 0
     
	 print "* Merging scores"
	 sys.stdout.write("Merge")
	 sys.stdout.flush()
     
	 maxpvalues = {} # sequence -> enzyme -> max pvalue
	 maxgtgscores = {} # sequence -> enzyme -> max gtg score
	 Fmult = {} # (sequence, enzyme) -> sum of scores KbKgKe
	 Fterms = {} 
	 Fmax = {} # (seq, enz) -> sum of max(Kb, Kg)Ke
	 Fmaxmax = {} # (seq, enz) -> max of max(Kb, Kg)Ke
	 # in keys saving all of the seq1 (query sequence id) 
	 keys = allseqs.keys()
	 keys.sort()
	 for seq1 in keys:
	     if seq1 not in Fmult:
		 Fmult[seq1] = {}
		 Fterms[seq1] = {}
		 maxpvalues[seq1] = {}
		 maxgtgscores[seq1] = {}
	     # save all seq2 for the same seq1
	     k2 = list(allseqs[seq1])
	     k2.sort()
	     for seq2 in k2:
		 if seq1 in blastpvalues and seq2 in blastpvalues[seq1]:
		     jointp = blastpvalues[seq1][seq2]
		 else:
		     jointp = "?"
		 # because no gtg right now, not concern
		 if seq1 in gtgscores and seq2 in gtgscores[seq1]:
		     gtgs = gtgscores[seq1][seq2]
		 else:
		     gtgs = "?"
		 #for each unknown input seq1, find all the ecs according seq2.
		 ecs = list(dbseqtoec[seq2])
		 ecstr = ",".join(ecs)
		 #below is the original one, I cancel gtgs for no gtg feature now
		 o.write("%s\t%s\t%s\t%s\t%s\n" % (seq1, seq2, ecstr, jointp, gtgs))
		 #below is the newline without gtg
		 #o.write("%s\t%s\t%s\t%s\n" % (seq1, seq2, ecstr, jointp))
		 if jointp == "?" or jointp < BASELINE_JOINTP:
		      jointp = BASELINE_JOINTP
		 if gtgs == "?" or gtgs < BASELINE_GTGSCORE:
		      gtgs = BASELINE_GTGSCORE
     
		 for ec in ecs:
		      if ec not in Fmult[seq1]:
			  Fmult[seq1][ec] = 0
			  Fterms[seq1][ec] = 0
		      Fmult[seq1][ec] += float(jointp) * float(gtgs)
		      Fterms[seq1][ec] += 1
		      if ec not in maxpvalues[seq1] or jointp > maxpvalues[seq1][ec]:
			  maxpvalues[seq1][ec] = jointp
		      if ec not in maxgtgscores[seq1] or gtgs > maxgtgscores[seq1][ec]:
			  maxgtgscores[seq1][ec] = gtgs
     
	     for ec in Fmult[seq1]:
		  kde = Fmult[seq1][ec] / Fterms[seq1][ec]
		  o2.write("%s\t%s\t%s\t%s\t%f\n" % (seq1, ec, maxpvalues[seq1][ec], maxgtgscores[seq1][ec], kde))
	 sys.stdout.write("\n")
