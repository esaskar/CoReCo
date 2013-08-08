#!/usr/bin/python
#
# McGregor's maximum common subgraph isomorphism (MCS) algorithm
# backtracking search algorithm, node-based
#
# basic algorithm is almost the same as VF2:
# - an recursive tree search procedure is used
# - search produces a search tree, where each node represents a partial mcs, 
#   with updated set of borders for this partial mcs
# - at each node, we try to match "first" atom on lhs border with all atoms on rhs border
# - a feasibility function checks that atoms match and that mapped neighbors of (n,m) are same
# - a mcs is found when search tree cannot be extended (i.e. we are at a leaf node)
#
# to change to edge-based mcs's, extend nodes with bond pairs instead of atom pairs
#


import sys, optparse, time
sys.path.append("/home/fs/mqheinon/duuni/python/")
from kegg import *
import time


class MCS(object):
	def mcs(self):
		pass


class Rascal(MCS):
	def __init__(self, g1, g2):
		self.g1 = g1
		self.g2 = g2
		self.iters = 0
		self.results = []
		self.verbose = False
		self.starttime = None
		self.endtime = None
	
	
	def modular_product_graph(self):
		self.mpg = kegg.MoleculeGraph()
		for a1,a2 in [(a1,a2) for a1 in self.g1.atoms for a2 in self.g2.atoms if a1.symbol == a2.symbol]:
			self.mpg.AddAtom()
			self.mpg.atoms[-1].a1 = a1
			self.mpg.atoms[-1].a2 = a2
		
		for v in self.mpg.atoms:
			for u in self.mpg.atoms:
				if v.a1.IsNeighbor(u.a1) and v.a1.IsNeighbor(u.a2):
					self.mpg.AddBond(u,v)
				elif not v.a1.IsNeighbor(u.a1) and not v.a1.IsNeighbor(u.a2):
					self.mpg.AddBond(u,v)
	
	def mcs(self):
		
		self.modular_product_graph()
		
		self.match()
	
	def match(self):
		pass
	




class McGregor(MCS):
	def __init__(self, g1, g2):
		self.g1 = g1
		self.g2 = g2
		self.iters = 0
		self.mcssize = 0
		self.mcslist = []
		self.verbose = False
		self.stop = False
		self.starttime = None
		self.runtime = None
	
	def match(self, s):
		self.iters += 1
		
		if self.stop or (self.runtime and time.time() - self.starttime > self.runtime):
			return
		
#		if self.iters > 100:
#			return
		
		if len(s) > self.mcssize:
			self.mcslist = [s.map] # save newest result as sole
			self.mcssize = len(s)
			
			if self.mcssize == len(self.g1) or self.mcssize == len(self.g2):
				self.stop = True
			
			if self.verbose:
				print "optimal mcs size now", len(s)
				print s
		elif len(s) == self.mcssize:
			self.mcslist.append(s.map)
			if self.verbose:
				print "+ mcs"
		
		if self.verbose:
			print " "*len(s), "state", s #, "[", lhs, rhs, "]"
		
		for n1, n2 in self.candidatepairs(s):
			if self.verbose:
				print " "*len(s), "cp [%s] [%s]" % (n1,n2), #state
			if self.feasible(s, n1,n2):
				if self.verbose:
					print
				
				ns = State()
				ns.map = s.map.copy()
				ns.map[n1] = n2
				
				ns.lhsborder = s.lhsborder[:]
				try: 
					ns.lhsborder.remove(n1) 
				except: 
					pass
				for ne in n1.GetAtomNeighbors():
					if ne not in ns.lhsborder and ne not in ns.map.keys():
						ns.lhsborder.append(ne)
				
				ns.rhsborder = s.rhsborder[:]
				try: 
					ns.rhsborder.remove(n2) 
				except: 
					pass
				for ne in n2.GetAtomNeighbors():
					if ne not in ns.rhsborder and ne not in ns.map.values():
						ns.rhsborder.append(ne)
				
				ns.lhsborder.sort()
				ns.rhsborder.sort()
				
				# only continue if border still contains stuff
				self.match(ns)
		
	
	# here need to check all pairs
	def candidatepairs(self, s):
		if s.lhsborder and s.rhsborder:
#			lhs = s.lhsborder[0]
			for lhs in s.lhsborder:
				for rhs in s.rhsborder:
					yield lhs,rhs
		
		# empty lhsborder!
		elif not s.lhsborder and not s.rhsborder:
			lhs = min( set(self.g1.atoms) - set(s.map.keys()) )
			for rhs in self.g2.atoms:
				if rhs not in s.map.values():
					yield lhs,rhs
	
	def feasible(self, s, n1, n2):
		if n1.symbol != n2.symbol:
			if self.verbose:
				print "ATOM MISMATCH"
			return False
		
		n_mapped_neighs = [ ne for ne in n1.GetAtomNeighbors() if ne in s.map.keys() ]
		m_mapped_neighs = [ ne for ne in n2.GetAtomNeighbors() if ne in s.map.values() ]
		
		if len(n_mapped_neighs) != len(m_mapped_neighs):
			if self.verbose:
				print "UNCONNECTED"
			return False
		
		for node in n_mapped_neighs:
			if not s.map[node] in m_mapped_neighs:
				if self.verbose:
					print "UNCONNECTED"
				return False
		
		return True

	def mcs(self):
		s = State() # empty state
		self.starttime = time.time()
		self.match(s)
		
		if self.runtime and time.time() - self.starttime > self.runtime:
			print "break ",
		
		return self.mcslist, self.mcssize


class State:
	def __init__(self):
		self.map = {}
		self.lhsborder = []
		self.rhsborder = []
	def __len__(self):
		return len(self.map)
	def __str__(self):
		s = "[" + ",".join([str(x.id) for x in self.map.keys()]) + "|" + ",".join([str(self.map[x].id) for x in self.map]) + "]"
		s += " [" + ",".join([str(x.id) for x in self.lhsborder]) + "]"
		s += " [" + ",".join([str(x.id) for x in self.rhsborder]) + "]" 
		return s

if __name__ == "__main__":
	
	parser = optparse.OptionParser()
	parser.set_usage("%prog [opts] mol1 mol2")
	parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose mode")
	
	(options, args) = parser.parse_args()
	
	if len(args) not in [1,2]:
		parser.error("give at least one molecule")
		parser.usage()
	
	if len(args) == 1:
		args.append(args[0])
	
	mol1 = Molecule(args[0])
	mol2 = Molecule(args[1])
	
#	mol1 = Molecule("C06161") # coclaurine, used also in EMCSS paper
#	mol2 = Molecule("C09084") # brucine, used also in EMCSS paper
	
	print "matching %s (size %d) and %s (size %d)" % (mol1.mol_id, len(mol1.graph), mol2.mol_id, len(mol2.graph))
	
	begintime = time.time()
	
	mcs = McGregor(mol1.graph, mol2.graph)
	mcs.verbose = options.verbose
	mcslist, mcssize = mcs.mcs()
	
	endtime = time.time()
	
	print "maximum mcs", mcssize
	print len(mcslist), "mcs's found"
	print "took %.2f secs and %d iters" % (endtime-begintime, mcs.iters)



