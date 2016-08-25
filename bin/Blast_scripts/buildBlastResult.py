#!/usr/bin/env python
"""
buildBlastResult.py - combine two-way blast results with EC numbers from UniProt
"""

import sys, datetime
# f2 is Organism->db;  f3 is db->Organism

def combineBlasts(f, f2, f3, o):
    print "Reading uniprot -> EC..."
    up2ec = {}
    for s in f:
        if s.startswith("#"):
            continue
        apu = s.strip().split() 
        uid = apu[0]
        pext = apu[1]
        if len(apu) == 3:
            ec = apu[2]
        else:
            ec = "?"
        up2ec[uid] = [pext, []]
        up2ec[uid][1] = ec.split(",")        
    blastres = {}

    # org -> uniprot blast
    print "Reading org -> uniprot blast results..."
    for s in f2:
        sseqid, qseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = s.strip().split("\t")        
        if sseqid not in blastres:
            blastres[sseqid] = {}
        blastres[sseqid][qseqid] = [evalue, float(bitscore), "?", "?"]

    # uniprot -> org blast
    print "Reading uniprot -> org blast results..."
    for s in f3:
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = s.strip().split("\t")
        # Discard result if sequence pair is not a hit in both directions
        if sseqid in blastres and qseqid in blastres[sseqid]:
            blastres[sseqid][qseqid][2] = evalue
            blastres[sseqid][qseqid][3] = float(bitscore)
       
    print "Writing joint results..."
    o.write("#Output of %s in %s\n" % (sys.argv[0], datetime.datetime.now()))
    o.write("# BlastFwd: target vs UniProt blast\n")
    o.write("# BlastRev: UniProt vs target blast\n")
    o.write("#TargetSeq UniProtSeq PExt ECS BlastFwdEValue BlastFwdScore BlastRevEValue BlastRevScore\n")
    keys1 = blastres.keys()
    keys1.sort()
    for s in keys1:
        keys2 = blastres[s].keys()

        def seqsort(x, y):
            xev1, xs1, xev2, xsc2 = blastres[s][x]
            yev1, ys1, yev2, ysc2 = blastres[s][y]
            return cmp(ys1, xs1)

        keys2.sort(seqsort)
        for q in keys2:
            if blastres[s][q][0] == "?" or blastres[s][q][2] == "?":
                continue
            upacc = q.split("|")[1]
            if upacc in up2ec:
                pext, ecs = up2ec[upacc]
            else:
                pext = "?"
                ecs = ["?"]
            for ec in ecs:
                o.write("%s\t%s\t%s\t%s\t%s\n" % (s, upacc, pext, ec, "\t".join(map(str, blastres[s][q]))))
   
if __name__ == "__main__":
    f = open(sys.argv[1])   # UniProt ID -> EC list
    f2 = open(sys.argv[2])  # blast res org -> uniprot
    f3 = open(sys.argv[3])  # blast res uniprot -> org
    o = open(sys.argv[4], "w")   # output
    combineBlasts(f, f2, f3, o)

