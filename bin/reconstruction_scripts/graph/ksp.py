#!/usr/bin/env python
#
# ksp.py - k shortest (simple!) paths
#
# Find k shortest paths using an external program
# that implements Yen's algorithm
#
# Time complexity is O(kn(m + n log n)) where
# k is the number of shortest paths enumerated,
# m and n are the number of edges and nodes.
#

import sys, os, subprocess
from graph import Graph

#KSP = '/group/home/icomic/software/ksp/ksp'
KSP = '~/programming/ksp/ksp'

namemap = {}
namemap2 = {}

def nnum():
    c = 0
    while 1:
        yield c
        c += 1

nn = nnum()

def mapname(s):
    if s not in namemap:
        nextix = nn.next()
        namemap[s] = nextix
        namemap2[nextix] = s
        nextix += 1
    return namemap[s]

def getindex(s):
    if s not in namemap:
        raise RuntimeError, "%s not in namemap" % (s)
    return namemap[s]

def getname(ix):
    if ix not in namemap2:
        raise RuntimeError, "%s not in namemap2" % (ix)
    return namemap2[ix]

def checkKSP():
    p = subprocess.Popen(KSP, shell = True, 
                         stdin=subprocess.PIPE, 
                         stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE, 
                         close_fds=True)

    return p.stdout.read().startswith("Usage")

def writeKSPinput(g, of):
    s = ""
    maxix = 0
    edges = 0
    for u in g.E:
        ui = mapname(u)
        if ui > maxix:
            maxix = ui
        for v in g.E[u]:
            vi = mapname(v)
            if vi > maxix:
                maxix = vi
            s += "%s\t%s\t%s\n" % (ui, vi, g.E[u][v])
            edges += 1
    
    of.write("%d\n\n%s" % (maxix + 1, s))

    #status("Wrote atom trace graph with %d nodes and %d edges" % (maxix, edges))

def findKSP(fn, numPaths, source, target):
    src = getindex(source)
    tgt = getindex(target)

    #print KSP, numPaths, fn, src, tgt
      
    #fin, fout, ferr = os.popen3("%s %d %s %s %s" % (KSP, numPaths, fn, src, tgt))

    p = subprocess.Popen("%s %d %s %s %s" % (KSP, numPaths, fn, src, tgt),
                         shell = True, 
                         stdin=subprocess.PIPE, 
                         stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE, 
              close_fds=True)

    #print fout.read()
    #print ferr.read()

    err = p.stderr.read()
    if err.startswith("The path"):
        print "Path %s -> %s does not exist" % (source, target)
        return None
    elif err.find("not found") != -1:
        print "Can't find %s!" % (KSP)
        return None

    #print namemap
    #print namemap2

    paths = []
    for s in p.stdout:
        vals = s.strip().split()
        id = vals[0]
        path = vals[1:]
        res = []
        for p in path:
            res.append(getname(int(p)))
        paths.append(res)

    return paths

def findKShortestPaths(g, k, s, t, tempfn = None):
    "Find k shortest paths from s to t in the graph g."

    if tempfn == None:
        ofn = os.tempnam()
    else:
        ofn = tempfn
    writeKSPinput(g, open(ofn, "w"))
    paths = findKSP(ofn, k, s, t)
    os.remove(ofn)
    return paths

def test():
    from randomGraph import randomUniformGraph
    g = randomUniformGraph(20, 1.0)
    paths = findKShortestPaths(g, 10, "N0", "N1")
    for path in paths:
        print path


if __name__ == "__main__":
    test()
