#!/usr/bin/env python

import sys, os, subprocess, math, random, re

import common

#MAX_BLAST_SCORE = 1-(1e-10)
ec_pa = re.compile("[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+")
LAPLACIAN_FRACTION = 0.20  # amount of noise added

def plot(outfn, orgname, blastmin = 0, blastmax = common.max_blast_joint_score(), format = "png"):

    if format == "png":
        formats = ",width=1024,height=1024,pointsize=24"
    else:
        formats = ""

    plotfn = "%s.R" % (outfn) 
    o = open(plotfn, "w")
    o.write("""
    bpos = read.table("%s.pos.blast")
    bneg = read.table("%s.neg.blast")
    gpos = read.table("%s.pos.gtg")
    gneg = read.table("%s.neg.gtg")

    blastpos = density(bpos[,1], from=%s, to=%s)
    blastneg = density(bneg[,1], from=%s, to=%s)
    gtgpos = density(gpos[,1], from=0, to=1)
    gtgneg = density(gneg[,1], from=0, to=1)

    name = c("BlastPos", "BlastNeg", "GTGPos", "GTGNeg")
    bw = c(blastpos[3], blastneg[3], gtgpos[3], gtgneg[3])

    write.table(cbind(blastpos$x, blastpos$y, blastpos$y / sum(blastpos$y)), file = "%s.blastpos")
    write.table(cbind(blastneg$x, blastneg$y, blastneg$y / sum(blastneg$y)), file = "%s.blastneg")
    write.table(cbind(gtgpos$x, gtgpos$y, gtgpos$y/sum(gtgpos$y)), file = "%s.gtgpos")
    write.table(cbind(gtgneg$x, gtgneg$y, gtgneg$y/sum(gtgneg$y)), file = "%s.gtgneg")
    write.csv(cbind(name, bw), file = "%s.bw")

    %s(file = "%s-all.%s"%s)

    par(mfrow=c(2,2), oma=c(2, 0, 3, 0), omi=c(0, 0, 0.8, 0))

    blastpos$y = blastpos$y / sum(blastpos$y)
    blastneg$y = blastneg$y / sum(blastneg$y)
    gtgpos$y = gtgpos$y / sum(gtgpos$y)
    gtgneg$y = gtgneg$y / sum(gtgneg$y)

    plot(blastpos, main = "P(BLAST | r = 1)")
    plot(blastneg, main = "P(BLAST | r = 0)")
    plot(gtgpos, main = "P(GTG | r = 1)")
    plot(gtgneg, main = "P(GTG | r = 0)")

    dev.off()
    """ % (outfn, outfn, outfn, outfn, 
           blastmin, blastmax, blastmin, blastmax, 
           outfn, outfn, outfn, outfn, outfn,
           format, outfn, format, formats))
    o.close()
    subprocess.call("R CMD BATCH %s" % (plotfn), shell = True)
    #subprocess.call("R --no-save <%s" % (plotfn), shell = True)

def add_laplacian_correction(nbpos, nbneg, ngpos, ngneg, obpos, obneg, ogpos, ogneg):
#def add_laplacian_correction(nbpos, nbneg, obpos, obneg):
    nbposn = int(LAPLACIAN_FRACTION * nbpos)
    nbnegn = int(LAPLACIAN_FRACTION * nbneg)
    ngposn = int(LAPLACIAN_FRACTION * nbpos)
    ngnegn = int(LAPLACIAN_FRACTION * nbneg)
    #print nposn, nnegn

    # add points uniformly to both score ranges
    for i in range(nbposn):
        bnoise = (1.0 * i / (nbposn - 1)) * common.max_blast_joint_score()
        obpos.write("%s\n" % (bnoise))
    for i in range(nbnegn):
        bnoise = (1.0 * i / (nbposn - 1)) * common.max_blast_joint_score()
        obneg.write("%s\n" % (bnoise))
    for i in range(ngposn):
        gnoise = (1.0 * i / (ngposn - 1))
        ogpos.write("%s\n" % (gnoise))
    for i in range(ngnegn):
        gnoise = (1.0 * i / (ngnegn - 1))
        ogneg.write("%s\n" % (gnoise))

def estimate(ecscoresfn, modelecsfn, outfn, use_raw_scores = False):

    f1 = open(ecscoresfn)   # .ecscores/BLAST, GTG
    f2 = open(modelecsfn)   # list of EC numbers in model

    orgname = modelecsfn[-4:]

    modelecs = set()
    for s in f2:
	if s.startswith("#"):
		continue
        modelecs.add(s.strip())
    #print modelecs
    obpos = open("%s.pos.blast" % (outfn), "w")
    obneg = open("%s.neg.blast" % (outfn), "w")
    ogpos = open("%s.pos.gtg" % (outfn), "w")
    ogneg = open("%s.neg.gtg" % (outfn), "w")
    header = "#MaxScore\n"
    obpos.write(header)
    obneg.write(header)
    ogpos.write(header)
    ogneg.write(header)
    nbpos = nbneg = ngpos = ngneg = 0
    #nbpos = nbneg = 0
    for line, s in enumerate(f1):
	#print s.split("\t")
        if s.startswith("#"):
            continue

        if use_raw_scores:
	    #EC BlastScore BlastQuerySeq BlastMatchSeq (.rawscore)
            seq1, seq2, ecs, maxblast, maxgtg = s.strip().split("\t")
	    #ecs, jointp, seq1, seq2 = s.strip().split("\t")
	    #print seq1, seq2, ecs, jointp
        else:
	    #.ecscore
	    #EC BlastScore BlastQuerySeq BlastMatchSeq GTGScore GTGQuerySeq GTGMatchSeq
            #print s
            ecs, maxblast, blastseq1, blastseq2, maxgtg, gtgseq1, gtgseq2 = s.strip().split("\t")
        if ecs == "?" or not ec_pa.match(ecs):
            continue
        
        if ecs in modelecs:
            if maxblast != "?":
                obpos.write("%s\n" % (maxblast))
                nbpos += 1
            if maxgtg != "?":
                ogpos.write("%s\n" % (maxgtg))
                ngpos += 1
	    #print ecs
        else:
            if maxblast != "?":
                obneg.write("%s\n" % (maxblast))
                nbneg += 1
            if maxgtg != "?":
                ogneg.write("%s\n" % (maxgtg))
                ngneg += 1

    #print "%d pos, %d neg" % (npos, nneg)

    add_laplacian_correction(nbpos, nbneg, ngpos, ngneg, obpos, obneg, ogpos, ogneg)
    #add_laplacian_correction(nbpos, nbneg, obpos, obneg)

    obpos.close()
    obneg.close()
    ogpos.close()
    ogneg.close()

    plot(outfn, orgname, 0, common.max_blast_joint_score(), "pdf")
    plot(outfn, orgname, 0, common.max_blast_joint_score(), "png")


if __name__ == "__main__":
    if len(sys.argv) > 4 and sys.argv[4] == "raw":
        use_raw_scores = True
    else:
        use_raw_scores = False
    estimate(sys.argv[1], sys.argv[2], sys.argv[3], use_raw_scores) #ecscoresfn, modelecsfn, outfn
