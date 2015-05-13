#!/usr/bin/env python
"""
Bayesian network
- Supports only discrete nodes
- Exact inference in polytrees with Pearl's algorithm

TODO:
- support for diagnostic evidence nodes with arbitrary distributions

Author: Esa Pitkanen (esa.pitkanen@helsinki.fi)
"""

import sys, math
import random
import traceback

from digraph import Digraph

VERBOSE = 0

class Message:

    LAMBDA_DONE = 0
    PI_DONE = 1
    LAMBDA = 2
    PI = 3

    def __init__(self, from_node, to_node, msg_type, msg):
        self.from_node = from_node
        self.to_node = to_node
        self.msg_type = msg_type
        self.msg = msg

def log2(x):
    return math.log(x) / math.log(2)

def logadd(logx, logy):
    """Logarithm of sum of two numbers given in log.

https://facwiki.cs.byu.edu/nlp/index.php/Log_Domain_Computations"""

    if logy > logx:
        logx, logy = logy, logx

    if logx == -float("inf"):
        return logx

    diff = logy - logx
    if diff < -53:  # does not make a difference at least in python 2.7.6
        return logx

    return logx + log2(1.0 + 2**diff)

class BayesianNetwork:

    EPS = 1e-300

    def __init__(self):
        self.pds = {}
        self.domain = {}
        # node id -> parent node id -> parent rank, e.g. {"v1" -> 0, "v2" -> 1}
        self.parent_order = {}

    def set_dag(self, dag):
        self.dag = dag

    def set_pd(self, node, parent_order, pd):
        self.domain[node] = pd.keys()
        self.parent_order[node] = {}
        for i in range(len(parent_order)):
            self.parent_order[node][parent_order[i]] = i
        log_pds = {}
        for x in pd:
            if type(pd[x]) == float:
                log_pds[x] = log2(pd[x])
            else:
                log_pds[x] = {}
                for y in pd[x]:
                    log_pds[x][y] = log2(pd[x][y])
        self.pds[node] = log_pds

    def compute_posterior(self, evidence):
        self.__la = {} 
        self.__pi = {}
        self.__la_msgs = {}
        self.__pi_msgs = {}
        self.__la_not_rec = {}
        self.__pi_not_rec = {}
        self.__la_not_sent = {}
        self.__pi_not_sent = {}

        log_evidence = {}

        for x in evidence:
            if type(evidence[x]) == float:
                log_evidence[x] = log2(max(BayesianNetwork.EPS, evidence[x]))
            else:
                log_evidence[x] = {}
                for y in evidence[x]:
                    log_evidence[x][y] = log2(max(BayesianNetwork.EPS, evidence[x][y]))

        messages = []
        for u in self.dag.parents:
            self.__la_msgs[u] = []
            self.__pi_msgs[u] = []
            self.__la_not_rec[u] = set(self.dag.children[u])
            self.__pi_not_rec[u] = set(self.dag.parents[u])
            self.__la_not_sent[u] = set(self.dag.parents[u])
            self.__pi_not_sent[u] = set(self.dag.children[u])
            if u in log_evidence:
                continue
            n = len(self.dag.parents[u])
            if n == 0:
               self.__pi[u] = self.pds[u]  # P(u = 1) prior
               messages.append(Message(None, u, Message.PI_DONE, None))
        for u in self.dag.children:
            if u in log_evidence:
                continue
            n = len(self.dag.children[u])
            if n == 0:
                self.__la[u] = {}
                for x in self.domain[u]:
                    self.__la[u][x] = 0 # log(1)
                messages.append(Message(None, u, Message.LAMBDA_DONE, None))

        for e in log_evidence:
            # evidence[e] = P(e = 1)
            self.__la[e] = log_evidence[e]
            self.__pi[e] = log_evidence[e]
            messages.append(Message(None, e, Message.LAMBDA_DONE, None))
            messages.append(Message(None, e, Message.PI_DONE, None))

        while len(messages) > 0:
            msg = messages.pop(0)
            if msg.msg_type == Message.LAMBDA_DONE:
                messages.extend(self.__handle_lambda_done(msg))
            elif msg.msg_type == Message.PI_DONE:
                messages.extend(self.__handle_pi_done(msg))
            elif msg.msg_type == Message.LAMBDA:
                messages.extend(self.__handle_lambda(msg))
            elif msg.msg_type == Message.PI:
                messages.extend(self.__handle_pi(msg))
            else:
                raise Exception("Unknown message %s" % (msg.msg_type))

        posterior = {}
        for u in self.dag.parents:
            posterior[u] = {}
            total = log2(BayesianNetwork.EPS)
            for x in self.domain[u]:
                total = logadd(total, self.__la[u][x] + self.__pi[u][x])
            for x in self.domain[u]:
                nbel = self.__la[u][x] + self.__pi[u][x] - total
                if VERBOSE:
                    print "BEL(%s = %s) = %s" % (u, x, nbel)
                posterior[u][x] = 2**nbel
        return posterior

    def __str__(self):
        ps = ""
        for u in self.dag.parents:
            if len(self.dag.parents[u]) == 0:
                for val in self.pds[u]:
                    ps += "\nP(%s=%s) = %s" % (u, val, self.pds[u][val])
            else:
                for val in self.pds[u]:
                    for pv in self.pds[u][val]:
                        ps += "\nP(%s=%s|pa(x)=%s) = %s" % (u, val, pv, self.pds[u][val][pv])
        s = "Bayesian network\n%s%s" % (self.dag, ps)
        return s

    def __handle_lambda(self, msg):
        u, v = msg.from_node, msg.to_node
        if VERBOSE:
            print "(l: %s -> %s) __handle_lambda = %s" % (u, v, msg.msg)
        if u not in self.__la_not_rec[v]:
            raise Exception("Duplication lambda msg %s -> %s" % (u, v))
        self.__la_not_rec[v].remove(u)
        self.__la_msgs[v].append(msg)

        msgs = []
        msgs.extend(self.__calculate_lambda(v))
        msgs.extend(self.__send_pi(v))
        return msgs

    def __calculate_lambda(self, u):
        if VERBOSE:
            print "\t_calculate_lambda", u
        if u not in self.__la and len(self.__la_not_rec[u]) == 0:
            self.__la[u] = {}  # value of u -> lambda(u)(value)
            for x in self.domain[u]:
                # e.g., 0 or 1
                la = 0 # log2(1.0)
                for msg in self.__la_msgs[u]:
                    la = la + msg.msg[x]
                self.__la[u][x] = la
            if VERBOSE:
                print "\tlambda(%s) = %s" % (u, self.__la[u])
            return [Message(None, u, Message.LAMBDA_DONE, None)]
        else:
            return []        

    def __get_pi_msg(self, u, v):
        pimsg = {}
        for x in self.domain[u]:
            val = self.__pi[u][x]
            for msg in self.__la_msgs[u]:
                if msg.from_node == v:
                    continue
                val = val + msg.msg[x]
            pimsg[x] = val
        return pimsg

    def __get_lambda_msg(self, Y, X):
        """Compute the lambda message Y -> X."""
        lamsg = {}
        for x in self.domain[X]:
            val = log2(BayesianNetwork.EPS)
            for y in self.domain[Y]:
                val2 = log2(BayesianNetwork.EPS)
                for parent_values in self.pds[Y][y]:
                    # e.g., parent_values == "0,1,0" or "true,false"
                    pv = parent_values.split(",")
                    ix = self.parent_order[Y][X]
                    if pv[ix] != x:
                        continue
                    val3 = self.pds[Y][y][parent_values]
                    for msg in self.__pi_msgs[Y]:
                        if msg.from_node == X:
                            continue
                        # jx = value of the parent 
                        jx = pv[self.parent_order[Y][msg.from_node]]
                        val3 = val3 + msg.msg[jx]
                    val2 = logadd(val2, val3)
                val2 += self.__la[Y][y]
                val = logadd(val, val2)
            lamsg[x] = val
        return lamsg                    

    def __send_pi(self, u):
        if VERBOSE:
            print "\t__send_pi", u, u in self.__pi, list(self.__la_not_rec[u]), list(self.__pi_not_sent[u])
        if u in self.__pi:
            if len(self.__la_not_rec[u]) == 1:
                v = list(self.__la_not_rec[u])[0]
                if v in self.__pi_not_sent[u]:
                    value = self.__get_pi_msg(u, v)
                    self.__pi_not_sent[u].remove(v)
                    return [Message(u, v, Message.PI, value)]
                else:
                    return []
            elif len(self.__la_not_rec[u]) == 0:
                msgs = []
                for v in self.__pi_not_sent[u]:
                    value = self.__get_pi_msg(u, v)
                    msgs.append(Message(u, v, Message.PI, value))
                self.__pi_not_sent[u].clear()
                return msgs
            else:
                return []
        else:
            return []

    def __calculate_pi(self, u):
        if VERBOSE:
            print "\t__calculate_pi", u
        if u not in self.__pi and len(self.__pi_not_rec[u]) == 0:
            self.__pi[u] = {}
            for x in self.domain[u]:
                # eg, 0 or 1
                total = log2(BayesianNetwork.EPS)
                for parent_values in self.pds[u][x]:
                    # e.g., "010" or "11"
                    pv = parent_values.split(",")
                    cp = self.pds[u][x][parent_values]
                    for msg in self.__pi_msgs[u]:
                        jx = pv[self.parent_order[u][msg.from_node]]
                        cp = cp + msg.msg[jx]
                    total = logadd(total, cp)
                self.__pi[u][x] = total
            if VERBOSE:
                print "\tpi(%s) = %s" % (u, self.__pi[u])
            return [Message(None, u, Message.PI_DONE, None)]
        else:
            return []

    def __send_lambda(self, u):
        if VERBOSE:
            print "\t__send_lambda", u
        if u in self.__la:
            if len(self.__pi_not_rec[u]) == 1:
                v = list(self.__pi_not_rec[u])[0]
                if v in self.__la_not_sent[u]:
                    value = self.__get_lambda_msg(u, v)
                    self.__la_not_sent[u].remove(v)
                    return [Message(u, v, Message.LAMBDA, value)]
                else:
                    return []
            elif len(self.__pi_not_rec[u]) == 0:
                msgs = []
                for v in self.__la_not_sent[u]:
                    value = self.__get_lambda_msg(u, v)
                    msgs.append(Message(u, v, Message.LAMBDA, value))
                self.__la_not_sent[u].clear()
                return msgs
        return []       

    def __handle_pi(self, msg):
        u, v = msg.from_node, msg.to_node
        if VERBOSE:
            print "(p, %s -> %s) __handle_pi: %s" % (u, v, msg.msg)
        if u not in self.__pi_not_rec[v]:
            raise Exception("Duplicate pi msg %s -> %s" % (u, v))
        self.__pi_not_rec[v].remove(u)
        self.__pi_msgs[v].append(msg)
        msgs = []
        msgs.extend(self.__calculate_pi(v))
        msgs.extend(self.__send_lambda(v))
        return msgs

    def __handle_lambda_done(self, msg):
        u = msg.to_node
        if VERBOSE:
            print "\t__handle_lambda_done", u, len(self.__pi_not_rec[u])
        return self.__send_lambda(u)

    def __handle_pi_done(self, msg):
        u = msg.to_node
        if VERBOSE:
            print "\t__handle_pi_done", u, len(self.__la_not_rec[u])
        return self.__send_pi(u)

def test2():
    dag = Digraph()
    dag.add_edge("v1", "v2")
    dag.add_edge("v2", "v3")
    bn = BayesianNetwork()
    bn.set_dag(dag)
    bn.set_pd("v1", [[0.2, 0.8]])
    bn.set_pd("v2", [[0.9, 0.1], [0.6, 0.4]])
    bn.set_pd("v3", [[0.5, 0.5], [0.3, 0.7]])

    evidence = {"v3" : {0 : 0.0, 1 : 1.0}}

    bn.compute_posterior(evidence)

def test():
    dag = Digraph()
    dag.add_edge("v1", "v2")
    dag.add_edge("v2", "v3")
    dag.add_edge("v2", "v4")

    bn = BayesianNetwork()
    bn.set_dag(dag)
    bn.set_pd("v1", [], {0 : 0.6, 1 : 0.4})
    bn.set_pd("v2", ["v1"], {0 : {"0" : 0.9, "1" : 0.1}, 1 : {"0" : 0.7, "1" : 0.3}})
    bn.set_pd("v3", ["v2"], {0 : {"0" : 0.1, "1" : 0.9}, 1 : {"0" : 0.8, "1" : 0.2}}) 
    bn.set_pd("v4", ["v2"], {0 : {"0" : 0.5, "1" : 0.5}, 1 : {"0" : 0.25, "1" : 0.75}})

    print bn

    evidence = {}
    bn.compute_posterior(evidence)

def test5():
    dag = Digraph()
    dag.add_edge("v1", "v2")
    dag.add_edge("v2", "v3")
    dag.add_edge("v2", "v4")

    bn = BayesianNetwork()
    bn.set_dag(dag)
    bn.set_pd("v1", [], {0 : 0.5, 1 : 0.5})
    bn.set_pd("v2", ["v1"], {0 : {"0" : 0.5, "1" : 0.5}, 1 : {"0" : 0.5, "1" : 0.5}})
    bn.set_pd("v3", ["v2"], {0 : {"0" : 0.5, "1" : 0.5}, 1 : {"0" : 0.5, "1" : 0.5}})
    bn.set_pd("v4", ["v2"], {0 : {"0" : 0.5, "1" : 0.5}, 1 : {"0" : 0.5, "1" : 0.5}})
    print bn

    evidence = {"v3" : {0 : 0.0, 1 : 1.0}, "v4" : {0 : 0.9, 1 : 0.1}}
    bn.compute_posterior(evidence)

def test3():
    dag = Digraph()
    dag.add_edge("v2", "v1")
    dag.add_edge("v3", "v1")

    bn = BayesianNetwork()
    bn.set_dag(dag)
    bn.set_pd("v1", ["v2", "v3"], {0 : {"11" : 0.8, "01" : 0.9, "10" : 0.5, "00" : 0.6}, 0 : {"11" : 0.2, "01" : 0.1, "10" : 0.5, "00" : 0.4}}) 
    bn.set_pd("v2", [], {1 : 0.1, 0 : 0.9})
    bn.set_pd("v3", [], {1 : 0.4, 0 : 0.6})

    print bn

    evidence = {}
    print bn.compute_posterior(evidence)

def test4():
    "v1 -> v2 -> v3"
    dag = Digraph()
    dag.add_edge("v1", "v2")
    dag.add_edge("v2", "v3")
    bn = BayesianNetwork()
    bn.set_dag(dag)
    bn.set_pd("v1", [], {"1" : 0.8, "0" : 0.2})
    bn.set_pd("v2", ["v1"], {"1" : {"1" : 0.4, "0" : 0.9}, "0" : {"1" : 0.6, "0" : 0.1}})
    bn.set_pd("v3", ["v2"], {"1" : {"1" : 0.7, "0" : 0.5}, "0" : {"1" : 0.3, "0" : 0.5}})
    evidence = {"v3" : {"1" : 1, "0" : 0}}
    print bn.compute_posterior(evidence)

def test6():
    dag = Digraph()
    dag.add_edge("v1", "v2")
    dag.add_edge("v2", "v3")
    bn = BayesianNetwork()
    bn.set_dag(dag)
    bn.set_pd("v1", [], {"v1_true" : 0.8, "v1_false" : 0.2})
    bn.set_pd("v2", ["v1"], {"v2_true" : {"v1_true" : 0.4, "v1_false" : 0.9}, "v2_false" : {"v1_true" : 0.6, "v1_false" : 0.1}})
    bn.set_pd("v3", ["v2"], {"v3_true" : {"v2_true" : 0.7, "v2_false" : 0.5}, "v3_false" : {"v2_true" : 0.3, "v2_false" : 0.5}})
    evidence = {"v3" : {"v3_true" : 1, "v3_false" : 0}}
    post = bn.compute_posterior(evidence)
    print post

def test7():
    dag = Digraph()
    dag.add_edge("v1", "v2")
    dag.add_edge("v2", "v3")
    bn = BayesianNetwork()
    bn.set_dag(dag)
    bn.set_pd("v1", [], {"v1_true" : 1e-300, "v1_false" : 1.0 - 1e-300})
    bn.set_pd("v2", ["v1"], {"v2_true" : {"v1_true" : 1e-300, "v1_false" : 1.0 - 1e-300}, "v2_false" : {"v1_true" : 1e-300, "v1_false" : 1.0 - 1e-300}})
    bn.set_pd("v3", ["v2"], {"v3_true" : {"v2_true" : 0.7, "v2_false" : 0.5}, "v3_false" : {"v2_true" : 0.3, "v2_false" : 0.5}})
    evidence = {"v3" : {"v3_true" : 1, "v3_false" : 0}}
    post = bn.compute_posterior(evidence)
    print "Posterior:", post

def testlog():
    x = 0.9
    y = 1e-320
    lx = log2(x)
    ly = log2(y)
    lz = logadd(lx, ly)
    print lx, ly, lz, 2**lz, "<==>", x, y, x + y

if __name__ == "__main__":
    test7()
