"""
Directed graph
"""

class Digraph(object):
    def __init__(self):
        self.children = {}
        self.parents = {}

    def add_edge(self, from_node, to_node):
        if from_node not in self.children:
            self.children[from_node] = set()
        if from_node not in self.parents:
            self.parents[from_node] = set()
        if to_node not in self.children:
            self.children[to_node] = set()
        if to_node not in self.parents:
            self.parents[to_node] = set()

        self.children[from_node].add(to_node)
        self.parents[to_node].add(from_node)

    def __str__(self):
        s = "Nodes: "
        k = self.children.keys()
        k.sort()
        for u in k:
            s += "%s, " % (u)
        s = s.rstrip(", ") + "\nEdges:\n"
        for u in k:
            k2 = list(self.children[u])
            k2.sort()
            s += "%s: " % (u)
            for v in k2:
                s += "%s, " % (v)
            s = s.rstrip(", ") + "\n"
        return s.rstrip("\n")
