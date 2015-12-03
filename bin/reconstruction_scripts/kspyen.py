#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
"""
ReTrace - find branching pathways in metabolic networks
Copyright (C) 2009 Esa Pitkänen
See file COPYING for license notice.

Implementation of Yen's k shortest simple (loopless) paths.

Running time is O(kn(m + n log n)) in a graph with n nodes
and m edges.

"""

#try:
#    import psyco # try to enable psyco JIT
#    psyco.full()
#except:
#    pass

import sys, os, random, time, math

from graph import randomGraph, graphGenerators
from eppstein import dijkstra, priodict

class Path(list):
    def __hash__(self):
        h = ""
        for e in self.__iter__():
            h += "%s*" % (e)
        return hash(h)

def getPath(s, t, P, G):
    cost = 0.0
    n = t
    p = Path()

    if len(P) == 0 or t not in P:
        return p, 0

    while n != s:
        cost += G.E[P[n]][n]
        p.append(n)
        n = P[n]
    p.append(n)
    p.reverse()
    return p, cost

def kspSimpleYen(G, K, s, t):
    assert(K > 0)
    X = priodict.priorityDictionary()  # paths
    removals = {}
    
    found = set()
    pathways = []

    (D, P) = dijkstra.Dijkstra(G.E, s, t)
    p, cost = getPath(s, t, P, G)
    p.cost = cost
    p.dix = 0
    X[p] = cost
    removals[p] = []
    found.add(p)

    k = 0
    for p in X:
        R = removals.pop(p)

        pathways.append(p)

        k += 1
        if k == K:
            break

        # Add to removals every edge from nodes v_{l-1}, l = |p|
        # on the path p, starting from the end of the path

        for (u, v, c) in R:
            del G.E[u][v]

        R2 = []
        i = p.dix
        while i < len(p) - 1:

            # Delete the next edge on the path
            pi = p[i]
            pii = p[i + 1]
            R2.append((pi, pii, G.E[pi][pii]))
            del G.E[pi][pii]

            (D, P) = dijkstra.Dijkstra(G.E, pi, t)
            p2, cost = getPath(pi, t, P, G)

            if len(p2) > 0:
                p3 = Path(p[0:i])
                p3.extend(p2)

                if p3 not in found:
                    p3.cost = p.cost + cost
                    p3.dix = i
                    X[p3] = p3.cost
                    removals[p3] = list(R2) 
                    found.add(p3)

            R3 = []
            for u in G.E[pi]:
                if u != pii:
                    R3.append((pi, u, G.E[pi][u]))

            for (u, v, c) in R3:
                del G.E[u][v]
            R2.extend(R3)
            i += 1

        # Restore graph

        for (u, v, c) in R2:
            G.E[u][v] = c

        for (u, v, c) in R:
            G.E[u][v] = c

    return pathways

def testKSP():
    G = graphGenerators.martinsPascoalExample()
    s = "N1"
    t = "N5"
    P = kspSimpleYen(G, 100, s, t)
    for p in P:
        print p

def meansd(L):
    mean = 1.0 * reduce(lambda x, y: x + y, L) / len(L)
    sum = 0.0
    for x in L:
        sum += (x - mean) * (x - mean)
    sd = math.sqrt(sum / len(L))
    return mean, sd


def testKSP2():
    edgetonoderatio = 1.5

    krange = range(50, 55, 5)
    nrange = range(200000, 1000000, 10000)
    rounds = 1
    rrange = range(rounds)

    o = sys.stdout

    s = 1
    t = 2

    for k in krange:
        for n in nrange:

            o.write("%d\t%d\t" % (k, rounds))
            o.flush()
            m = int(edgetonoderatio * n)
            G = graphGenerators.randomSparseGraph(n, m)
            o.write("%d\t%d\t" % (G.numNodes(), G.numEdges()))
            o.flush()
            times = []

            for r in rrange:
                V = list(G.V)
                s = random.choice(V)
                t = random.choice(V)

                st = time.time()

                paths = kspSimpleYen(G, k, s, t)

                et = time.time()

                times.append(et - st)

            mean, sd = meansd(times)
            o.write("%f\t%f\n" % (mean, sd))
            o.flush()

def removeCyclicPaths(P):
    Q = []
    for p in P:
        nodes = set()
        cyclic = False
        for n in p:
            if n in nodes:
                cyclic = True
                break
            else:
                nodes.add(n)
        if not cyclic:
            Q.append(p)
    return Q

if __name__ == "__main__":
    testKSP2()
   
