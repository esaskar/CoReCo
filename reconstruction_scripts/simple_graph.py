import scc # Tarjan's algorithm

def is_cyclic(E):
    sccs = scc.strongly_connected_components(E)
    for comp in sccs:
        if len(comp) > 1:
            return True
    return False

def find_reachable(E, srcs):
    Q = list(srcs)
    visited = set(srcs)
    while len(Q) > 0:
        q = Q.pop(0)
        if q in E:
            for u in E[q]:
                if u not in visited:
                    Q.append(u)
                    visited.add(u)
    return visited

def backtrack(Erev, nodes):
    Q = list(nodes)
    visited = set(nodes)
    ends = set()
    edges = {}
    #print ">> NODES", nodes, " EDGES", Erev
    while len(Q) > 0:
        q = Q.pop(0)
        if q in Erev and len(Erev[q]) > 0:
            for u in Erev[q]:
                if u not in edges:
                    edges[u] = set()
                if q not in edges:
                    edges[q] = set()
                edges[u].add(q)
                if u not in visited:
                    Q.append(u)
                    visited.add(u)
        else:
            ends.add(q)
    return ends, edges

def scc_compress(E):
    sccs = scc.strongly_connected_components(E)

    #print "*** SCCS", sccs

    node_to_comp = {}
    comp_to_node = {}
    C = {}
    Crev = {}
    for i, comp in enumerate(sccs):
        C[i] = {}
        Crev[i] = {}
        comp_to_node[i] = set()
        for u in comp:
            node_to_comp[u] = i
            comp_to_node[i].add(u)

    for u in E:
        cu = node_to_comp[u]
        for v in E[u]:
            cv = node_to_comp[v]
            if cu != cv:
                C[cu][cv] = 1
                Crev[cv][cu] = 1

    return C, Crev, node_to_comp, comp_to_node

ccc = 0

def scc_backtrack(Erev, nodes):
    Crev, C, node_to_comp, comp_to_node = scc_compress(Erev)

    #print "scc_backtrack: NODES", nodes
    #print "scc_backtrack: Crev", Crev
    #print "scc_backtrack: ntc", node_to_comp
    
    start_comps = set()
    unmapped_nodes = set()
    for u in nodes:
        #print u, "->", node_to_comp[u]
        if u in node_to_comp:
            start_comps.add(node_to_comp[u])
        else:
            unmapped_nodes.add(u)

    end_comps, used_edges = backtrack(Crev, start_comps)
    
    end_nodes = set(unmapped_nodes)
    for c in end_comps:
        # arbitrarily choose one node from each compartment
        cnodes = set(comp_to_node[c])
        #print "COMP", c, "->", cnodes

        end_nodes.add(cnodes.pop())

    #print "END_NODES", end_nodes

    #global ccc
    #ccc += 1
    #if ccc == 50:
    #    import sys
    #    sys.exit()

    return end_nodes

def scc_backtrack_sources(Erev, nodes, sources):
    Crev, C, node_to_comp, comp_to_node = scc_compress(Erev)

    start_comps = set()
    end_nodes = set()
    for u in nodes:
        if u in node_to_comp:
            start_comps.add(node_to_comp[u])
        else:
            end_nodes.add(u)

    end_comps, used_edges = backtrack(Crev, start_comps)
    
    for c in end_comps:
        # arbitrarily choose one node from each compartment that does not contain
        # a source
        has_source = False
        for u in comp_to_node[c]:
            if u in sources:
                has_source = True
                break
        if not has_source:
            cnodes = set(comp_to_node[c])
            end_nodes.add(cnodes.pop())

    return end_nodes

def quickdraw(E, fn):
    import subprocess, tempfile
    f = tempfile.NamedTemporaryFile()
    f.write("digraph bla {\n")
    for u in E:
        for v in E[u]:
            f.write("   %s->%s;\n" % (u, v))
    f.write("}\n")
    f.flush()
    subprocess.call("dot -Tsvg %s -o %s" % (f.name, fn), shell = True)
    f.close()
