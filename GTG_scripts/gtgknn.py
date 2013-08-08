#!/usr/bin/env python
"""
gtgknn.py - Retrieve k nearest GTG neighbors for a set of sequences with GTG features.

Input: 
- GTG features for sequences from gtg_attributes_Jun09.py
- GTG -> sequence index
- k - number of nearest hits 
Output:
- k nearest hits for each input sequence

Author: Esa Pitkanen (esa.pitkanen@cs.helsinki.fi)
"""

import sys, array

from buildGTGindex import ADDR_SIZE

# Index is inconveniently large, use files instead
# def readIndex(f):
#     gtg2seq = {}
#     for s in f:
#         gtg, seqs = s.strip().split("\t")
#         gtg2seq[int(gtg)] = map(int, seqs.split(","))
#     return gtg2seq

def findMatches(seqgtgs, ptrf, dataf):
    #print seqgtgs
    matchseqs = {}
    for gtg in seqgtgs:
        addr = (gtg - 1) * ADDR_SIZE

        try:
            ptrf.seek(addr)
        except:
            raise
        ptra = array.array("I")
        assert(ptra.itemsize == ADDR_SIZE)
        try:
            ptra.fromfile(ptrf, 1)  # pointers to this and next gtg data
            #print "POS", ptrf.tell(), "read:", ptra[0]

        except:
            raise
        daddr = ptra[0]

        ptra = array.array("I")
        try:
            ptra.fromfile(ptrf, 1)
            nextaddr = ptra[0]
            #print "POS", ptrf.tell(), "read:", ptra[0]
        except:
            # attempted to access last gtg record
            ptrf.seek(0, 2) # eof
            nextaddr = ptrf.tell()
            
        nseqs = (nextaddr - daddr) / ADDR_SIZE
        #print gtg, addr, daddr, nextaddr, nseqs
        assert(nseqs >= 0)

        if nseqs > 0:
            try:
                dataf.seek(daddr)
            except:
                raise
            seqa = array.array("I")
            seqa.fromfile(dataf, nseqs)
            for seq in seqa:
                if seq not in matchseqs:
                    matchseqs[seq] = 1
                else:
                    matchseqs[seq] += 1
    return matchseqs
            
def reportKNearestMatches(seqid, seqgtgcount, matchseqs, k, seq2up, of):
    #print "***\t%s\t%d" % (seqid, seqgtgcount)

    if len(matchseqs) == 0:
        of.write("%s\t-\t-\t-\t-\n" % (seqid))
        return

    mseqs = matchseqs.keys()
    mseqs.sort(lambda x, y: cmp(matchseqs[y], matchseqs[x]))
    for seq in mseqs[:k]:
        #print "%s\t%d\t%s\t%s\t%s\t%.3f" % (seqid, seqgtgcount, seq, seq2up[seq], matchseqs[seq], 1.0 * matchseqs[seq] / seqgtgcount)
        if seq in seq2up:
            upid = seq2up[seq]
        else:
            upid = "?"
        of.write("%s\t%d\t%s\t%s\t%s\t%.3f\n" % (seqid, seqgtgcount, upid, seq, matchseqs[seq], 1.0 * matchseqs[seq] / seqgtgcount))

def readSeq2Up(f):
    seq2up = {}
    for s in f:
        if s.startswith("#"):
            continue
        seq, up = s.strip().split("\t")
        seq2up[int(seq)] = up
    return seq2up

def gtgknn(f, ixfn, seq2upf, k, of):
    #print "Reading gtg->seq index..."
    #gtg2seq = readIndex(ixf)
    
    ptrf = open("%s.ptr" % (ixfn), "rb")
    dataf = open("%s.data" % (ixfn), "rb")
    assert(ptrf)
    assert(dataf)

    print "Reading seq->uniprot index..."
    seq2up = readSeq2Up(seq2upf)
    currentseq = None
    seqgtgs = set()
    print "Finding k-NN..."
    of.write("#SeqId SeqGTGs MatchSeq MatchNID MatchGTGs MatchFrac\n")
    for s in f:
        seqid, pos, aa, seqref, pos2, aaix, gtg = s.strip().split("\t")
        gtg = int(gtg)
        if currentseq == None:
            currentseq = seqid
        if currentseq != seqid:
            #print currentseq, len(seqgtgs)
            matchseqs = findMatches(seqgtgs, ptrf, dataf)
            reportKNearestMatches(currentseq, len(seqgtgs), matchseqs, k, seq2up, of)
            currentseq = seqid
            seqgtgs = set()
        seqgtgs.add(gtg)
    matchseqs = findMatches(seqgtgs, ptrf, dataf)
    reportKNearestMatches(seqid, len(seqgtgs), matchseqs, k, seq2up, of)

if __name__ == "__main__":
    f = open(sys.argv[1])   # input gtg features
    #ixf = open(sys.argv[2]) # gtg -> seq index
    ixfn = sys.argv[2]  # gtg -> seq index, "gtg2index.CAAx"
    seq2up = open(sys.argv[3]) # seq -> uniprot accession index "nids.up"
    k = int(sys.argv[4])
    of = open(sys.argv[5], "w")
    gtgknn(f, ixfn, seq2up, k, of)

