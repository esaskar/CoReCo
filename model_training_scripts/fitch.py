
import sys, random

def parsimony(parents, children, models, all_ecs):
    # init
    nodes = {}
    Q = []
    child_not_done = {}
    root = None
    for u in parents:
        nodes[u] = {}
        if parents[u] == u:
            root = u
        if u not in children:
            Q.append(u)
        for ec in all_ecs:
            if u in models:
                if ec in models[u]:
                    nodes[u][ec] = 1
                else:
                    nodes[u][ec] = 0
            else:
                child_not_done[u] = len(children[u])

    # bottom up
    while len(Q) > 0:
        q = Q.pop(0)
        if q in children:
            nodes[q] = {}
            lc, rc = children[q][0], children[q][1]
            for ec in all_ecs:
                lv, rv = nodes[lc][ec], nodes[rc][ec]
                if lv == 0:
                    if rv == 1:
                        x = 2
                    else:
                        x = 0
                elif lv == 1:
                    if rv == 0:
                        x = 2
                    else:
                        x = 1
                else:
                    x = rv
                nodes[q][ec] = x
        if q in parents:
            p = parents[q]
            child_not_done[p] -= 1
            if child_not_done[p] == 0:
                Q.append(p)

    # top down
    Q = list(children[root])
    for ec in all_ecs:
        v = nodes[root][ec]
        if v == 2:
            nodes[root][ec] = random.randint(0, 1)

    num_add = 1
    num_del = 1
    num_no_add = 1
    num_no_del = 1

    mutations = {}
    while len(Q) > 0:
        q = Q.pop(0)
        parent = parents[q]
        ladd = ldel = lnadd = lndel = 1
        for ec in all_ecs:
            v = nodes[q][ec]
            pv = nodes[parent][ec]
            if v == 2:
                v = nodes[q][ec] = pv
            if v == 2:
                raise Exception("bad value")
            if v == 1:
                if pv == 1:
                    num_no_del += 1
                    lndel += 1
                else:
                    num_add += 1
                    ladd += 1
            else:
                if pv == 1:
                    num_del += 1
                    ldel += 1
                else:
                    num_no_add += 1
                    lnadd += 1
        key = "%s-%s" % (parent, q)
        mutations[key] = [ladd, ldel, lnadd, lndel] 
        #print q, parent, 1.0 * ladd / (ladd + lnadd), 1.0 * ldel / (ldel + lndel), ladd, ldel, lnadd, lndel
        if q in children:
            Q.extend(children[q])

    #print num_add, num_del, num_no_add, num_no_del
    return num_add, num_del, num_no_add, num_no_del, mutations
