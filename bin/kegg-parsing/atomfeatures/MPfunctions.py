#!/usr/bin/python
# -*- coding: iso-8859-15 -*-
#
# Message passing algorithm uses 5 different functions to control the
# algorithm
#
# - node_initfunc() = the initial value for nodes
#                     e.g. for atom distribution, set [null null C null ...]
# - node_feat()     = the g() value for nodes
#                     e.g. its own atom, rest zeros
# - msg_feat()      = the value for messages, used only at the beginning
#                     e.g. messages are initialized to source-node value
# - node_sumfunc()  = summarization function for nodes,
#                     e.g. sum of incoming messages plus own node_feat
# - msg_sumfunc()   = summarization function for messages,
#                     e.g. min-function for minimum distance computation
#


import Message


#####################################################################
#
# initialization functions for summarised node values
#
#
#def node_initfunc(graph, node, node_feat):
#	return node_feat(graph, node)
#
#def emptyset_node_initfunc(graph, node, node_feat):
#	return []
#
#def zero_node_initfunc(graph, node, node_feat):
#	return 0
#
#def bond_distribution_initfunc(graph, node, node_feat):
#	return Message.Vector(len(graph.atoms), "", node.id-1, "", node)
#
#
#class MorganInitfunc:
#	def __call__(self, graph, node, node_feat):
#		return Message.Vector(len(graph.atoms), None, node.id-1, None, node)
#
#
#
#####################################################################
#
# node features
#

class ExtdataNodeFeat:
	def __init__(self, extdata):
		self.extdata = extdata
	
	def __call__(self, graph, node):
		return Message.Vector(len(graph.atoms), None, node.id-1, self.extdata[node.id].vec, node)
	
class DistanceNodeFeat:
	def __call__(self, graph, node):
		return Message.Vector(len(graph.atoms), None, node.id-1, 0)
	

class AtomNodeFeat:
	def __call__(self, graph, node):
		return Message.Vector(len(graph.atoms), "", node.id-1, node.symbol)



def node_feat(graph, node):
	# apply the feature value to 'node_value',
	#  e.g. set own value to 1 and rest 0 ( [0 ... 1 ... 0] )
	return Message.Vector(len(graph.atoms), 0, node.id-1, 1, node)


def distance_node_feat(graph, node):
	return Message.Vector(len(graph.atoms), None, node.id-1, 0, node)

def atom_node_feat(graph, node):
	return Message.Vector(len(graph.atoms), "", node.id-1, node.symbol, node)

def carbon_node_feat(graph, node):
	return Message.Vector(len(graph.atoms), 99, node.id-1, (0 if node.symbol == "C" else 99))

def bond_node_feat(graph, node):
	# return bond types
	
	bonds = node.GetBondNeighbors()
	nf = Message.Vector(len(graph.atoms), "", 0, "", node)
	
	for b in bonds:
		if b.source == node:
			nf.max(Message.Vector(len(graph.atoms), "", b.target.id-1, str(b.type), b.target))
		elif b.target == node:
			nf.max(Message.Vector(len(graph.atoms), "", b.source.id-1, str(b.type), b.source))
	
	return nf

def extdata_node_feat(graph, node):
	# import 'fs' dictionary from FeatureGenerator
	from FeatureGenerator import ext
#	print ext, ext.data
	# use that as the unique value for this node
	return Message.Vector(len(graph.atoms), None, node.id-1, fs[node.id].vec, node)



class MorganNodeFeat:
	def __call__(self, graph, node):
		return Message.Real(len(node.bonds), node)
#		return Message.Real(1)


#####################################################################
#
# message features
#


def msg_feat(graph, edge, behind_node, node_feat):
	# apply this for initial values for msgs
	return node_feat(graph, behind_node)

def bond_msg_feat(graph, edge, behind_node, node_feat):
#	return node_feat(graph, behind_node)
	return Message.Vector(len(graph.atoms), "", behind_node.id-1, str(edge.type), behind_node)

class MorganMsgFeat:
	def __call__(self, graph, edge, behind_node, node_feat):
#		return node_feat(graph, behind_node)
		return Message.Real(len(behind_node.bonds), behind_node)



#####################################################################
#
# summarization functions for nodes
#


def node_sumfunc(graph, node, node_feat, node_val, msgs, d):
	val = node_feat
	for m in msgs:
		val += m
	return val


def min_distance_node_sumfunc(graph, node, node_feat, node_val, msgs, d):
	val = node_feat
	for m in msgs:
		val.min(m + 1)
	return val

def ring_node_sumfunc(graph, node, node_feat, node_val, msgs, ringsize):
	nid = node.id-1
	dists = [m.vec[nid] for m in msgs]
	
	if not dists:
		return False
	
	largest = max(dists)
	
	if largest:
		largest += 1
	
	if largest == ringsize:
		return True
	return False
	


def atom_distribution_node_sumfunc(graph, node, node_feat, node_val, msgs, d):
	val = node_feat
	for m in msgs:
		val.max(m)
	
	s = "(" + ''.join(sorted(val.vec)) + "".join( [ "," + str(m.node.id) + "=" + "".join(sorted(m.vec)) for m in msgs ] ) + ")"
	
	return s
#	return val

def atom_distribution_exclude_node_sumfunc(graph, node, node_feat, node_val, msgs, d):
	val = Message.Vector(len(graph.atoms)).fill("")
	for m in msgs:
		val.max(m)
	
#	print node_feat
#	print val
	s = "(" + ''.join(sorted(val.vec)) + "".join( [ "," + str(m.node.id) + "=" + "".join(sorted(m.vec)) for m in msgs ] ) + ")"
	
	return s



def bond_distribution_node_sumfunc(graph, node, node_feat, node_val, msgs, d):
	val = node_feat
	
	for m in msgs:
		val.max(m)
	
	mainres = "".join(sorted([elem for elem in val.vec if elem]))
	dirres = "".join([",%s=%s" % (str(m.node.id), "".join(sorted([elem for elem in m.vec if elem]))) for m in msgs ])
	
	s = "(" + mainres + dirres + ")"
	
	return s


def wiener_node_sumfunc(graph, node, node_feat, node_val, msgs, d):
	gsize = len(graph.atoms)
	reached = [[] for x in range(gsize)]
	sumvalues = [0 for x in range(gsize)]
	
	# compute subfunctions for each message
	for m in msgs:
		m_id = m.node.id-1
		# append node_value to m
		m.max(node_feat)
		m = m.vec
		
		# compute reached = indices of non-null values
		reached[m_id] = [index for index,value in enumerate(m) if value is not None]
		
		# for each dist-list, take those which are reached in that direction
		for i in reached[m_id]:
			v = [value for index,value in enumerate(m[i]) if index in reached[m_id]]
			sumvalues[m_id] += sum(v)
	
	# compute whole function for the node
	val = node_feat
	
#	print "original value", val
	vid = val.node.id-1
	for m in msgs:
#		print "value to be added", m
		val.max(m)
#		print "value after adding", val
#	print "resulting", val
	
	reached[vid] = [index for index,value in enumerate(val.vec) if value is not None]
	
#	print "reached[vid]", reached[vid]
	
#	print "XXX", reached[vid]
	for i in reached[vid]:
#		print val.vec[i]
		v = [value for index,value in enumerate(val.vec[i]) if index in reached[vid]]
#		print i+1, val.vec[i], v, sum(v)
		sumvalues[vid] += sum(v)
	
	# divide by two
	sumvalues = [x/2 for x in sumvalues]
	
	values_str = ["," + str(ind+1) + "=" + str(val) for ind,val in enumerate(sumvalues) if val != 0 and ind != vid]
	
	s = "(" + str(sumvalues[vid]) + "".join(values_str) + ")"
	
	return s
#	return sumvalues[vid]
	
	


def wiener2_node_sumfunc(graph, node, node_feat, node_val, msgs, d):
	# produce a min-message over all
	res = "("
#	print "wiener"
	if node.id == 11:
		print "before", msgs[0]
	val = node_feat
#	print "start", val
	for m in msgs:
#		print m
		val.min(m + 1)
#	print val
	res += str(val.sum()) + ","
	
#	if node.id == 11:
#		print "whole sum", val.sum()
#		print "before", msgs[0]
	
	
	for e in node.GetBondNeighbors():
		val = None
#		if node.id == 11:
#			print "source", e.source.id, "target", e.target.id
		if e.source != node:
			val = e.msgforward[-2]
#			if node.id == 11:
#				print "edge as it is", val
			# take only smallest distances
			for m in msgs:
				val.min_or_none(m)
#			if node.id == 11:
#				print "after min_or_none", val
			val += 1
#			if node.id == 11:
#				print "after +1", val
			res += str(e.source.id) + ":" + str(val.sum())
			
#			if node.id == 11:
#				print val.sum()
		else:
			val = e.msgbackward[-2]
#			if node.id == 11:
#				print "edge as it is", val
			
			# take only smallest distances
			for m in msgs:
				val.min_or_none(m)
#			if node.id == 11:
#				print "after min_or_none", val
#			print val
			val += 1
#			print val
#			if node.id == 11:
#				print "after +1", val
			res += str(e.target.id) + ":" + str(val.sum())
#			if node.id == 11:
#				print val.sum()
		res += ", "
	node_feat(graph,e.target),
	res = res.strip(", ")
	res += ")"
#	print "end", res
	return res
	
#	print "end", val
#	print "sum", val.sum()
	# take a sum over it
	
	return val.sum()

#	return val.sum()


class MorganNodeSumfunc:
	def __call__(self, graph, node, node_feat, node_val, msgs, d):
		if d == 0:
			return "(0)"
#			return "(" + str(node_feat) + ")"
		
		val = Message.Real(0)
		for m in msgs:
			val += m

		s = "(" + str(val) + "".join( [ "," + str(m.node.id) + "=" + str(m.value) for m in msgs ] ) + ")"
		
		return s
		

#####################################################################
#
# message summarization functions
#


def msg_sumfunc(graph, behind_node, behind_node_feat, msgs):
	# summarizes the messages and the behind node
	val = behind_node_feat
	for m in msgs:
		val += m
	return val

def min_distance_msg_sumfunc(graph, behind_node, behind_node_feat, msgs):
	val = behind_node_feat
	# msgs + 1
	for m in msgs:
		val.min(m + 1)
	return val

def max_distance_msg_sumfunc(graph, behind_node, behind_node_feat, msgs):
	val = behind_node_feat
	# msgs + 1
	for m in msgs:
		val.max(m + 1)
	return val

def atom_distribution_msg_sumfunc(graph, behind_node, behind_node_feat, msgs):
	val = behind_node_feat
	for m in msgs:
		val.max(m)
	return val

def bond_distribution_msg_sumfunc(graph, behind_node, behind_node_feat, msgs):
	# joka atomilla on naapurikaarien tyyppijoukko, esim "121" tai "111"
	# n�ist� tulee ottaa maximit?

	val = behind_node_feat
	for m in msgs:
		val.max(m)
	
	return val



class MorganMsgSumfunc:
	def __call__(self, graph, behind_node, behind_node_feat, msgs):
		val = Message.Real(0, behind_node)
		for m in msgs:
			val += m
		return val



#######################################



class MorganEndfunc:
	def __call__(self, nodes):
		if len(nodes[0].val) < 11:
			return False
		
		return True


def summarize(gsize, node, iter):
	node.val.append(Message.Vector(gsize, node))
	for z in node.GetAtomNeighbors():
		node.val[iter].update( z.val[iter-1], int.__add__ )

def summarize_rings(node, node_val, msgs):
	# from messages, count whether any of them has a distance "6" for itself
	for m in msgs:
		if m.vec[node.id-1] == 6:
			return True
	return False




