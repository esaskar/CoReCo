from graph import Graph
from random import random

def completeGraph(n, p):
    """p disregarded"""
    g = Graph()
    for i in range(n):
        g.addNode("N%d" % (i))
    for u in g.V:
        for v in g.V:
            g.addEdge(u, v, random())    
    return g

def smallExampleGraph():
    g = Graph()
    for i in range(1, 6):
        g.addNode("N%d" % (i))
    g.addEdge("N1", "N5", 0)
    g.addEdge("N2", "N1", 4)
    g.addEdge("N2", "N4", 1)
    g.addEdge("N3", "N2", 1)
    g.addEdge("N3", "N4", 2)
    g.addEdge("N4", "N1", 1)
    g.addEdge("N5", "N3", 1)
    g.addEdge("N5", "N2", 3)
    return g
