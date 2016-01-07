#!/usr/bin/env python

class MetabolicNetwork:

    def __init__(self, reactions):
        self.G = {}
        self.Grev = {}
        for i, rid in enumerate(reactions):
            r = reactions[rid]
            if rid not in self.G:
                self.G[rid] = {}
                self.Grev[rid] = {}
            for mol in r.substrates:
                if mol not in self.G:
                    self.G[mol] = {}
                self.G[mol][rid] = 1
                self.Grev[rid][mol] = 1
            for mol in r.products:
                self.G[rid][mol] = 1
                if mol not in self.Grev:
                    self.Grev[mol] = {}
                self.Grev[mol][rid] = 1

class AtomGraph:

    def __init__(self, reactions, cost_fun):
        """Build an atom graph. Argument reactions is a iterable collection of atommap.Reaction instances."""
        self.E = {}
        self.Erev = {}
        self.edge_reactions = {}
        for r in reactions:
            self.add_reaction_edges(reactions[r], cost_fun(r))

    def add_reaction_edges(self, re, cost):	
        for sa in re.maps:
            for ta in re.maps[sa]:
		 self.add_edge(sa, ta, cost, [re.id])
		

    def delete_reaction_edges(self, re):
        for sa in re.maps:
            for ta in re.maps[sa]:
                assert(re.id in self.edge_reactions[sa][ta])
		self.edge_reactions[sa][ta].remove(re.id)
                if len(self.edge_reactions[sa][ta]) == 0:
                     del self.E[sa][ta]
                     del self.Erev[ta][sa]
		


    def scale_edge_weights(self, amp):
        """Scale edge weights w' = w/x where x is the fraction of source atoms transferred to target."""
        S = {}
        srcs = {}
        #outedges = {}
        #counts = {}
        for u in self.E:
            mol1, atom1 = amp.get_mol_and_atom(u)
            if mol1 not in S:
                #counts[mol1] = {}
                S[mol1] = {}
                srcs[mol1] = {}
                #outedges[mol1] = 0
            for v in self.E[u]:
                mol2, atom2 = amp.get_mol_and_atom(v)
                if mol2 not in S[mol1]:
                    #counts[mol1][mol2] = 1
                    S[mol1][mol2] = set([(u, v)])
                    srcs[mol1][mol2] = set([u])
                else:
                    #counts[mol1][mol2] += 1
                    S[mol1][mol2].add((u, v))
                    srcs[mol1][mol2].add(u)
                #outedges[mol1] += 1

        for mol1 in S:
            for mol2 in S[mol1]:
                n = len(srcs[mol1][mol2])
                #x = 1.0 * n / outedges[mol1]
                x = 1.0 * n / len(amp.atom_types[mol1])
                assert(0 <= x <= 1)
                for u, v in S[mol1][mol2]:
                    self.E[u][v] = 1.0 * self.E[u][v] / x
                    self.Erev[v][u] = 1.0 * self.Erev[v][u] / x

    def add_edge(self, u, v, cost, edge_res):
        """Add edge (u, v) to graph and update cost if the new cost is lower than previous, or the edge is new."""
	if u not in self.E:
            self.E[u] = {}
            self.edge_reactions[u] = {}
	if v not in self.E[u] or self.E[u][v] > cost:
            self.E[u][v] = cost
	if v not in self.edge_reactions[u]:
            self.edge_reactions[u][v] = set()
        self.edge_reactions[u][v].update(edge_res)
	if v not in self.Erev:
            self.Erev[v] = {}
        if u not in self.Erev[v] or self.Erev[v][u] > cost:
            self.Erev[v][u] = cost

    def delete_edge(self, u, v):
        if u not in self.E or v not in self.E[u]:
            raise Exception("Invalid edge (%s,%s)" % (u, v))
        del self.E[u][v]
        del self.Erev[v][u]
        del self.edge_reactions[u][v]

    def delete_node(self, u):
        if u in self.E:
            for v in self.E[u]:
                del self.Erev[v][u]
            for v in self.Erev[u]:
                del self.E[v][u]
            del self.E[u]
            del self.Erev[u]
            del self.edge_reactions[u]

    def __str__(self):
        s = ""
        for u in self.E:
            for v in self.E[u]:
                s += "%s -> %s %s\n" % (u, v, list(self.E[u][v]))
        return s.rstrip("\n")

    def __repr__(self):
        return "AtomGraph (%d nodes)" % (len(self.E))
        
def test():
    ag = AtomGraph([])
    ag.add_edge("n1", "n0", ["r1"])
    ag.add_edge("n1", "n2", ["r1"])
    ag.add_edge("n2", "n3", ["r1"])
    ag.add_edge("n3", "n2", ["r1"])
    ag.add_edge("n3", "n4", ["r1"])
    ag.add_edge("n3", "n5", ["r1"])
    #ag.delete_edge("n2", "n3")
    ag.delete_node("n2")
    #print ag
    #print ag.Grev

if __name__ == "__main__":
    test()
