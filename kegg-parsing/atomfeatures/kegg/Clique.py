#!/usr/bin/python

import os, sys, time, optparse, random
import kegg




def bron_kerbosch(g):
	
	return bk(set(), set(g.nodes), set())
	

def bk(R, P, X):
#	print map(str, R), map(str,P), map(str,X)
	
	# P and X empty -> clique
	if len(P) == 0 and len(X) == 0:
#		print " maximal clique", len(R), map(str, R)
		return
	
	u = choose_pivot(P|X)
	
	for v in P:
		if v in u.neighs:
			continue
		
		bk( R|set([v]), P&set(v.neighs), X&set(v.neighs) )
		
		P.remove(v)
		X.add(v)



def choose_pivot(PX):
	return iter(PX).next()
	


NODES = 1000
EDGERATIO = 0.10

g = kegg.Graph()

for i in range(NODES):
	g.addnode()

for n1,n2 in [(n1,n2) for n1 in g.nodes for n2 in g.nodes if n1 != n2]:
	if random.random() < EDGERATIO:
		g.addedgeto(n1,n2)

print "graph done"

bron_kerbosch(g)









