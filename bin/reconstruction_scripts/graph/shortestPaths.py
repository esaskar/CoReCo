#!/usr/bin/env python
#
# shortestPaths.py
#

import sys
from graph import Graph
from eppstein.dijkstra import Dijkstra

def dijkstra(g, s, e = None):
    return Dijkstra(g.E, s, e)

def hashlist(l):
    s = ""
    items = list(l)
    items.sort()
    for (u, v, w) in items:
        s += "%s-%s*" % (u, v)
    return hash(s)

def hashpath(l):
    s = ""
    for i in l:
        s += "%s=" % (i)
    return hash(s)

def KSPDijkstra(g, k, s, e = None):
    G = Graph()

    for v in g.V:
        G.addNode(v)

    for u in g.E:
        for v in g.E[u]:
            G.addEdge(u, v, g.E[u][v])

    queue = [set()]
    paths = []

    visited = set()
    found = set()

    ix = 0
    c = 0
    while len(queue) > 0 and c < k:
        
        ix += 1
        removals = queue.pop(0)
        #print "Remove %d edges" % (len(removals))
        #print "Removing", list(removals)

        h = hashlist(removals)
        if h in visited:
            continue
        visited.add(h)

        for (u, v, w) in removals:
            #print "Removed: (%s,%s)" % (u, v)
            G.removeEdge(u, v)

        (D, P) = Dijkstra(G.E, s, e)

        for (u, v, w) in removals:
            G.addEdge(u, v, w)

        if e not in D:
            continue

        node = e
        path = [e]
        while node != s:
            node = P[node]
            path.insert(0, node)
        
        h = hashpath(path)
        if h in found:
            continue
        
        c += 1

        #print c, path
        paths.append(path)
        found.add(h)

        for i in range(len(path) - 1):
            newRemovals = set()
            newRemovals.update(removals)
            newRemovals.add((path[i], path[i + 1], G.E[path[i]][path[i + 1]]))
            
            if hashlist(newRemovals) in visited:
                continue

            #print "Queued:", list(newRemovals)

            queue.append(newRemovals)

    return paths

def testksp():
    from random import random
    g = Graph()
    n = 3 * 3
    for u in range(n):
        g.addNode(u)

    g.addEdge(0, 1, 1)
    g.addEdge(1, 2, 1)
    g.addEdge(0, 3, 1)
    g.addEdge(1, 4, 1)
    g.addEdge(2, 5, 1)
    g.addEdge(3, 4, 1)
    g.addEdge(4, 5, 1)
    g.addEdge(3, 6, 1)
    g.addEdge(4, 7, 1)
    g.addEdge(5, 8, 1)
    g.addEdge(6, 7, 1)
    g.addEdge(7, 8, 1)

    n = 10000
    g = Graph()
    for u in range(n):
        g.addNode(u)
    p = 0.3
    for u in range(n):
        for v in range(n):
            if random() < p:
                g.addEdge(u, v, random())

    print "Finding ksp..."
    sys.stdout.flush()
    KSPDijkstra(g, 10, 0, 8)
    

def test():
    from random import random
    g = Graph()
    n = 20
    for u in range(n):
        g.addNode(u)
    p = 0.1
    for u in range(n):
        for v in range(n):
            if random() < p:
                g.addEdge(u, v, random())

    (D, P) = Dijkstra(g.E, 0)
    print D
    print P

if __name__ == "__main__":
    #test()
    testksp()
