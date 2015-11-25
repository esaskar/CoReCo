#!/usr/bin/env python
"""
tree.py
"""

def components(A):
    Q = A.keys()
    C = {}
    c = 0
    while len(Q) > 0:
        u = Q.pop()
        if u not in C:
            C[u] = c
            c += 1
        for v in A[u]:
            if v not in C:
                Q.append(v)
                C[v] = C[u]
    return c, C

def treeToStr(A):
    s = ""
    for u in A:
        for v in A[u]:
            if u <= v:
                s += "%s\t%s\n" % (u, v)
    return s

def rootedTreeToStr(parent, children, root):
    s = ""
    Q = [root]
    while len(Q) > 0:
        q = Q.pop(0)
        if q != root:
            s += "%s\t%s\n" % (parent[q], q)
        if q in children:
            Q.extend(children[q])
    return s


def adjMapToEdgeList(A):
    E = []
    for u in A:
        for v in A[u]:
            if u <= v:
                E.append((u, v))
    return E

def edgeListToAdjMap(E):
    A = {}
    for u, v in E:
        if u not in A:
            A[u] = set()
        if v not in A:
            A[v] = set()
        A[u].add(v)
        A[v].add(u)
    return A

def findParentsAndChildren(A, root):
    parent = {}
    children = {}
    Q = []
    parent[root] = root
    for u in A[root]:
        if u != root:
            Q.append((root, u))
            parent[u] = root
    while len(Q) > 0:
        (u, v) = Q.pop(0)
        if u not in children:
            children[u] = []
        children[u].append(v)
        for w in A[v]:
            if w not in children:
                Q.append((v, w))
                parent[w] = v
    return parent, children

def findLeaves(A):
    """Return leaf nodes in tree. Specify tree as an adjacency map."""
    L = set()
    for u in A:
        if len(A[u]) == 1:
            L.add(u)
    return L

def numberOfInternalEdges(A):
    if components(A)[0] != 1:
        raise RuntimeError("Graph contains more than one component")
    return len(A) - len(findLeaves(A)) - 1

def deleteNonTerminalShoots(terminals, A):
    Q = findLeaves(A)
    #print "Before:", A
    #print "Leaves:", Q
    while len(Q) > 0:
        q = Q.pop()
        if q not in terminals:
            if len(A[q]) == 1:
                v = A[q].pop()
                del A[q]
                A[v].remove(q)
                Q.add(v)
            elif len(A[q]) == 0:
                del A[q]
    #print "After:", A
    return A

def contractNonJunctions(adj, protected = set()):
    A = adj.copy()
    Q = [findLeaves(A).pop()]
    V = set()
    while len(Q) > 0:
        q = Q.pop(0)
        V.add(q)
        for u in A[q]:
            if u not in V:
                Q.append(u)
        if len(A[q]) == 2 and q not in protected:
            u, v = list(A[q])
            A[u].add(v)
            A[v].add(u)
            A[q].remove(u)
            A[q].remove(v)
            A[u].remove(q)
            A[v].remove(q)
    return A

def expandShoots(adj, nodes, firstExtraIndex = -1):
    """Expand selected internal nodes to shoots."""
    A = {}

    mi = firstExtraIndex  # first index of new nodes
    for u in adj:
        A[u] = set()
        if u > mi:
            mi = u

    n = mi + 1
    replaced = {}
    for u in nodes:
        if len(adj[u]) > 1:  # is internal?
            # u becomes the shoot and is replaced with replaced[u]
            replaced[u] = n
            A[n] = set()
            n += 1

    adds = []
    for u in adj:
        ru = u
        if u in replaced:
            ru = replaced[u]
        for v in adj[u]:
            rv = v
            if v in replaced:
                rv = replaced[v]
            adds.append((ru, rv))
  
    for u, v in adds:
        A[u].add(v)
        A[v].add(u)
    
    for u in replaced:
        # shoot u is connected to replaced[u] only
        A[u] = set([replaced[u]])
        A[replaced[u]].add(u)

    return A, replaced

def getNewick(adj):
    for u in adj:
        if len(adj[u]) > 1:
            break
    assert(len(adj[u]) > 1)
    s = ""
    Q = [(u, False), (u, True)]
    visited = set([u])
    # dfs traversal of internal nodes
    while len(Q) > 0:
        (u, preorder) = Q.pop()
        if preorder:
            s += "("
            for v in adj[u]:
                if len(adj[v]) > 1:
                    if v not in visited:
                        Q.append((v, False))
                        Q.append((v, True))
                        visited.add(v)
                else:
                    s += "%s," % (str(v))
        else:
            s = s.rstrip(",") + "),"
    s = s.rstrip(",") + ";"
    return s

class NewickVisitor:
    def pre_visit_tree(self, x):
        pass
    def pre_visit_edge(self, t, b, l, n):
        node1 = self.mapnode(t)
        node2 = self.mapnode(n)
        if node1 not in self.A:
            self.A[node1] = set()
        if node2 not in self.A:
            self.A[node2] = set()
        self.A[node1].add(node2)
        self.A[node2].add(node1)
    def visit_leaf(self, l):
        pass
    def post_visit_edge(self, t, b, l, n):
        pass
    def post_visit_tree(self, x):
        pass
    def mapnode(self, n):
        n = str(n).strip("'")
        if not n.startswith("("):
            return n
        if n not in self.nodemap:
            self.nodemap[n] = self.nextnode
            self.nextnode += 1
        return "N%s" % (self.nodemap[n])
    def getAdjMap(self):
        return self.A
    def __init__(self):
        self.A = {}
        self.nodemap = {}
        self.nextnode = 1

def newickToTree(s):
    from newick import parse_tree #http://www.brics.dk/~mailund/newick.html
    t = parse_tree(s)
    v = NewickVisitor()
    t.dfs_traverse(v)
    A = v.getAdjMap()
    return A

def renameTree(A, names):
    B = {}
    for u in A:
        B[names[u]] = {}
        for v in A[u]:
            B[names[u]][names[v]] = A[u][v]
    return B

def testNewickToTree():
    print "=== testNewickToTree() ==="
    s = "(uma,((cne,dpch),(spo,((((ure,cim),(aor,(ani,afm))),(((mgr,ncr),fgr),ssl)),(yli,(((pgu,dha),(ctp,cal)),((cgr,sce),(kla,ago))))))));"
    s = "(uma,((cne,dpch),(spo,(((aor,(ani,afm)),((mgr,ncr),fgr)),(yli,((dha),((cgr,sce),(kla,ago))))))));"
    s = "(uma,((cne,dpch),(spo,(((aor,(ani,afm)),((mgr,ncr),fgr)),(yli,(dha,((cgr,sce),(kla,ago))))))));"

    A = newickToTree(s)
    E = adjMapToEdgeList(A)
    #print A
    #print E
    import taxa2image
    labels = {}
    for u in A:
        labels[u] = u
    taxa2image.graphToImage("test1.ps", E, labels)
    print "Check test1.ps"

def testNewick():
    E = [(1,3),(2,3),(3,4),(4,5),(4,6),(6,8),(6,7),(4,9),(9,10),(9,11),(6,12),(12,13),(12,14)]
    A = edgeListToAdjMap(E)
    s = getNewick(A)
    print s
    E = [(1,2),(4,2),(2,3),(3,5),(3,6)]
    A = edgeListToAdjMap(E)
    s = getNewick(A)
    print s

def testExpandShoots():
    import taxa2image
    labels = {}
    E = [(1,3),(2,3),(3,4),(4,5),(4,6),(6,8),(6,7)]
    A = edgeListToAdjMap(E)
    nodes = set([4,3])
    T = adjMapToEdgeList(expandShoots(A, nodes))
    for u in A:
        labels[u] = u
    print T
    taxa2image.graphToImage("test1.ps", T, labels)

    nodes = set([3,4,6])
    E = [(1,3),(2,3),(3,4),(4,5),(4,6),(6,8),(6,7),(4,9),(9,10),(9,11),(6,12),(12,13),(12,14)]
    E = [(1,3),(2,3),(3,4),(4,5),(4,6),(6,8),(6,7)]
    A = edgeListToAdjMap(E)
    T = adjMapToEdgeList(expandShoots(A, nodes))
    labels = {}
    for u in A:
        labels[u] = u
    print T
    taxa2image.graphToImage("test2.ps", T, labels)

    nodes = set([3,4,6,9,10,12])
    E = [(1,3),(2,3),(3,4),(4,5),(4,6),(6,8),(6,7),(4,9),(9,10),(9,11),(6,12),(12,13),(12,14)]
    A = edgeListToAdjMap(E)
    T = adjMapToEdgeList(expandShoots(A, nodes))
    labels = {}
    for u in A:
        labels[u] = u
    print T
    taxa2image.graphToImage("test3.ps", T, labels)

def testContract():
    import taxa2image
    E = [(0,7),(7,8),(8,9),(1,9),(9,11),(2,10),(3,10),(10,11),(11,14),(13,14),(13,6),(12,13),(4,12),(5,12)]
    A = edgeListToAdjMap(E)

    labels = {}
    for u in A: labels[u] = u
    taxa2image.graphToImage("test.ps", E, labels)

    T = contractNonJunctions(A)
    T = adjMapToEdgeList(T)
    labels = {}
    for u in A: labels[u] = u
    taxa2image.graphToImage("test2.ps", T, labels)

def testDeleteShoots():
    print "Testing deleteNonTerminalShoots()..."
    E = [(0,7),(7,8),(8,9),(1,9),(9,11),(2,10),(3,10),(10,11),(11,14),(13,14),(13,6),(12,13),(4,12),(5,12)]
    A = edgeListToAdjMap(E)
    terminals = set()
    assert(len(deleteNonTerminalShoots(terminals, A)) == 0)
    
    A = edgeListToAdjMap(E)
    terminals = set([9, 11, 14, 12, 6])
    assert(len(deleteNonTerminalShoots(terminals, A)) == 6)
    print "Success!"

def testComponents():
    print "Testing components()..."
    E = [(0,7),(7,8),(8,9),(1,9),(9,11),(2,10),(3,10),(10,11),(11,14),(13,14),(13,6),(12,13),(4,12),(5,12),(15,16),(16,17),(16,18),(19,20)]
    A = edgeListToAdjMap(E)
    c, C = components(A)
    assert(c == 3)
    E = [(0,7),(7,8),(8,9),(1,9),(9,11),(2,10),(3,10),(10,11),(11,14),(13,14),(13,6),(12,13),(4,12),(5,12)]
    A = edgeListToAdjMap(E)
    c, C = components(A)
    assert(c == 1)
    E = [(u, u + 1) for u in range(1, 10, 2)]
    A = edgeListToAdjMap(E)
    c, C = components(A)
    assert(c == 5)
    print "Success!"

def testIntEdges():
    E = [(1,3),(2,3),(3,4),(4,5),(4,6),(4,9),(9,8),(9,7),(5,10),(10,11)]
    A = edgeListToAdjMap(E)
    assert(numberOfInternalEdges(A) == 4)
    print "testIntEdges: Success!"

if __name__ == "__main__":
    #testDeleteShoots()
    #testComponents()
    #testContract()
    #testExpandShoots()
    #testNewick()
    #testIntEdges()
    testNewickToTree()
