
from graph import Graph
from random import random

def randomUniformGraph(n, p):
    g = Graph()
    for i in range(1, n + 1):
        g.addNode(i)
    for u in g.V:
        for v in g.V:
            if random() < p:
                g.addEdge(u, v, random())
    return g

def randomGraph(n, m, randomf = random):
    p = 1.0 * m / (n * n)

    g = Graph()
    for i in range(1, n + 1):
        g.addNode(i)
    for u in g.V:
        for v in g.V:
            if random() < p:
                g.addEdge(u, v, randomf())
    return g
