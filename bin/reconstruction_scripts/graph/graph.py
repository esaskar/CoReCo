#!/usr/bin/env python
#
# graph.py
#

class Graph:

    def __init__(self, V = [], E = {}):
        #self.V = list(V)
        self.V = set(V)
        self.E = dict(E)

    def numNodes(self):
        return len(self.V)

    def numEdges(self):
        return self.__countEdges()

    def addNode(self, u):
        #self.V.append(u)
        self.V.add(u)
        self.E[u] = {}

    def addEdge(self, u, v, w = 1):
        if u not in self.V:
            #self.V.append(u)
            self.V.add(u)
        if v not in self.V:
            #self.V.append(v)
            self.V.add(v)
        if u not in self.E:
            self.E[u] = {}
        if v not in self.E: 
            self.E[v] = {}
        self.E[u][v] = w
        
    def removeEdge(self, u, v):
        del self.E[u][v]

    def removeNode(self, u):
        for a in self.E:
            if u in self.E[a]:
                del self.E[a][u]
        if u in self.E:
            del self.E[u]
        self.V.remove(u)

    def removeNodes(self, U):
        for u in U:
            self.removeNode(u)

    def __str__(self):
        s = ""
        for u in self.E:
            for v in self.E[u]:
                s += "%s -> %s: %s\n" % (u, v, self.E[u][v])
        return s.strip()

    def __countEdges(self):
        c = 0
        for u in self.E:
            c += len(self.E[u])
        return c
                

def test():
    V = [1,2,3]
    E = {1:{2:1.0}, 2:{3:2.0}}
    g = Graph(V, E)
    print g
    print
    g.addNode(4)
    g.addEdge(1, 4, 3.0)
    print g
    print
    g.removeNode(1)
    print g
    print
    try:
        g.removeNode(1)
        print g
        print
    except KeyError, e:
        print "Ok"
    

if __name__ == "__main__":
    test()
