#!/usr/bin/env python
"""
buildGTGindex.py - Build a binary index for GTG -> sequence data

Author: Esa Pitkanen (esa.pitkanen@cs.helsinki.fi)
"""

import sys, array

from linebuffer import Linebuffer

ADDR_SIZE = 4
GTG_SIZE = 4

def buildIndex(f, ixfn):
    ptrf = open("%s.ptr" % (ixfn), "wb")
    dataf = open("%s.data" % (ixfn), "wb")
    lb = Linebuffer(f)
    nextgtg = 1
    nextaddress = 0
    while 1:
        if (nextgtg % 100000) == 0:
            print nextgtg
        s = lb.getLine()
        if s == "": #eof
            break
        caa, seqs = s.strip().split("\t")
        caa = int(caa)
        ptra = array.array("I") 
        assert(ptra.itemsize == ADDR_SIZE)
        if caa == nextgtg:
            seqs = map(int, seqs.split(","))
            ptra.append(nextaddress)
            nextaddress += len(seqs) * GTG_SIZE
            da = array.array("I")
            for seq in seqs:
                da.append(seq)
            da.tofile(dataf)                
            lb.advance()
        elif caa > nextgtg:
            ptra.append(nextaddress)   # skip this feature, not in data
        else:
            print nextgtg, s
            assert(0)
        ptra.tofile(ptrf)
        nextgtg += 1

if __name__ == "__main__":
    f = open(sys.argv[1])        # gtg2seq.txt
    ixfn = sys.argv[2]           # gtg2seq.index
    buildIndex(f, ixfn)
