#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import Message
import string, glob, optparse, sys, os, re
from MPfunctions import *


def AlgorithmA(graph, maxiter):
	"""AlgorithmA is the basic version of message passing algorithm
	where there are only factor nodes which send messages to eachother.
	
	Each node sends its information to all neighbors every round
	"""
	
	gsize = len(graph.atoms)
	
	# initialize all node values to [0 0 0 1 0 0]
	for n in graph.atoms:
		n.val = [ Message.Vector(gsize, n) ]
	
	# Loop for 'maxiter' times
	for d in range(1,maxiter):
		# In each loop,
		#  go through the nodes and summarize over neighbours +
		#  the original value of that node
		for n in graph.atoms:
			# g(x,z_1...z_n)
			
			n.val.append( Message.Vector(gsize,n) )
			# summarize f(z_1) ... f(z_n) to the g(x)
			for z in n.GetNodeNeighbors():
				n.val[d].update( z.val[d-1] )



def AlgorithmB(graph, node_feat, msg_feat, node_sumfunc, msg_sumfunc, limit_rounds = False, verbose = False, maxiter = 100, selfloops = False):
	gsize = len(graph.atoms)
	# default value
#	maxiter = 100
	
	if verbose:
		print "Starting MP algorithm"
		print "initializing node vectors"
	
	# initialize node values with 'node_feat' function
	#  (e.g. distance function: [99...99 0 99...99],
	#        where i'th is 0, rest infinity )
	for n in graph.atoms:
#		n.val = [ node_initfunc(graph, n, node_feat) ]
		n.val = [ node_sumfunc(graph, n, node_feat(graph, n), None, [], 0)]
#		n.val = [ node_feat(graph, n) ]
		
		if verbose:
			print " node", n.id, n.val[0]

	# initialize messages with 'msg_feat' function
	#  (e.g. msg = source-value )
	for e in graph.bonds:
		e.msgforward  = [ msg_feat(graph, e, e.source, node_feat) ]
		e.msgbackward = [ msg_feat(graph, e, e.target, node_feat) ]
#		e.msgforward  = [ node_feat(graph, e.source) ]
#		e.msgbackward = [ node_feat(graph, e.target) ]
	
		if verbose:
			print " edge", e.id, "messages", e.msgforward[0], e.msgbackward[0]
	
	
	# loop for 'maxiter' rounds
	for d in range(1, maxiter+1):
		
		if verbose:
			print "round", d
		
		messages_identical = True
		nodes_identical = True
		
		
		# optional ...
		for n in graph.atoms:
			# make room for this round
			n.val.append(Message.Vector(gsize))
			neighs = n.GetBondNeighbors()[:]
			
			# gather the messages
			msgs = [ (z.msgforward[d-1] if z.target == n else z.msgbackward[d-1]) for z in neighs ]
			
			# sum over all messages (and possibly the node itself)
			n.val[d] = node_sumfunc(graph, n, node_feat(graph,n), n.val[d-1], msgs, d)
			
			# check whether node values have changed since last iteration
			if n.val[d] != n.val[d-1]:
				nodes_identical = False
			
			if verbose:
				print " node", n.id, n.val[d]
		
		
		# append room for this round values
		for e in graph.bonds:
			e.msgforward.append(Message.Vector(gsize))
			e.msgbackward.append(Message.Vector(gsize))
		
		# compute new messages from neighbors
		for e in graph.bonds:
			
			# <<< forward cases >>>
			# sum over neighbors
			neighs = e.source.GetBondNeighbors()[:]
			if not selfloops:
				neighs.remove(e)
			msgs = [ (z.msgforward[d-1] if z.target == e.source else z.msgbackward[d-1]) for z in neighs ]
			
			# use sumfunc to go over the messages
			e.msgforward[d] = msg_sumfunc(graph, e.source, node_feat(graph,e.source), msgs)
			if e.msgforward[d] != e.msgforward[d-1]:
				messages_identical = False
			
			
			# <<< backward cases >>>
			# sum over neighbors
			neighs = e.target.GetBondNeighbors()[:]
			if not selfloops:
				neighs.remove(e)
			msgs = [ (z.msgforward[d-1] if z.target == e.target else z.msgbackward[d-1]) for z in neighs ]
			
			e.msgbackward[d] = msg_sumfunc(graph, e.target, node_feat(graph,e.target), msgs)
			if e.msgbackward[d] != e.msgbackward[d-1]:
				messages_identical = False
			
			if verbose:
				print " edge", e.id, "messages", e.msgforward[d], e.msgbackward[d]
		

		if limit_rounds:
			if messages_identical and nodes_identical:
#				print "executed", d, "rounds only"
				break
	
	# end



def WriteFeaturesList(mol, graph, feature_name, context, join=True):
	f = open("features.txt", "a")
	s = ""
	
	for c in context:
		for n in graph.atoms:
			data = n.val[c].vec
			data.sort()
			if join:
				s = mol + ":" + str(n.id-1) + ":" + feature_name + "_" + str(c) + "=" +''.join(data) + "\n"
			else:
				s = mol + ":" + str(n.id-1) + ":" + feature_name + "_" + str(c) + "=" + str(data) + "\n"
	
			f.write(s)
	
	f.close()


def WriteFeaturesValue(mol, graph, feature_name, context):
	f = open("features.txt", "a")
	s = ""
	
	for c in context:
		for n in graph.atoms:
			data = n.val[c]
			s = mol + ":" + str(n.id-1) + ":" + feature_name + "_" + str(c) + "=" + str(data) + "\n"
			
			f.write(s)
	
	f.close()


def WriteStats(graph):
	f = open("log.txt","w")
	
	f.write("Node values ...\n")
	for n in graph.atoms:
		f.write(str(n.id) + "\n")
		
		for i in range(len(n.val)):
			f.write( " " + str(i) + "\t" + str(n.val[i]) + "\n")

	try:
		f.write("\nEdge messages ...\n")
		for e in graph.bonds:
			f.write(str(e.id) + "\n")
			
			for i in range(len(e.msgfor)):
				f.write(" " + str(i) + "\t" + str(e.msgfor[i]) +" \t" + str(e.msgbac[i]) + "\n")
	except:
		pass
	
	f.close()




