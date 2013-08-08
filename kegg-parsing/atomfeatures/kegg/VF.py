#!/usr/bin/python
# -*- coding: iso-8859-15 -*-
#
# implements VF (sub)graph isomorphism algorithm
#

import sys, optparse, time
sys.path.append("/home/fs/mqheinon/duuni/python/")
from kegg import *

ISOMORPH = 0
SUBISOMORPH = 1
MAPPINGISOMORPH = 2

class GraphVF:
	ISOMORPH = 0
	SUBISOMORPH = 1
	MAPPINGISOMORPH = 2

	def __init__(self, graph1, graph2):
		self.graph1 = graph1
		self.graph2 = graph2
		self.itercount = 0
		self.ISOMODE = None
		self.solutions = []
		self.verbose = None
		self.isofound = False
		self.boolmode = False
	
	def candidatepairs(self, state):
		if state.lhsborder and state.rhsborder:
			rhs = state.rhsborder[0]
			for lhs in state.lhsborder:
				yield lhs,rhs
		
		# borders are empty!
		elif not state.lhsborder and not state.rhsborder:
			m = min( set(self.graph2.atoms) - set(state.map.values()) )
			for n in self.graph1.atoms:
				if n not in state.map.keys():
					yield n,m
	
	def feasible(self, state, n, m):
		# the regions of G1 and G2 are divided into three areas:
		# - (MR) mapped regions
		# - (BR) adjacent border of the mapped region
		# - (RR) remote regions
		#
		# semantic rule:
		#  atoms have to match
		#
		if n.symbol != m.symbol:
			if self.verbose:
				print "ATOM MISMATCH"
			return False
		
		# connectedness-rule
		#
		# see that M[neighbors(n)] == neighbors(m), i.e. that n's neighbors in M are m's neighbors
		#
		n_mapped_neighs = [ ne for ne in n.GetBondNeighbors() if ne.other(n) in state.map.keys() ]
		m_mapped_neighs = [ ne for ne in m.GetBondNeighbors() if ne.other(m) in state.map.values() ]
		
		if len(n_mapped_neighs) != len(m_mapped_neighs):
			if self.verbose:
				print "UNCONNECTED"
			return False
		
		for edge in n_mapped_neighs:
			src = edge.other(n)
			b = state.map[src].GetBond(m)
			if not b or b.type != edge.type:
				if self.verbose:
					print "UNCONNECTED"
				return False
		
		
		
#		for node in n_mapped_neighs:
#			if not state.map[node] in m_mapped_neighs:
#				if self.verbose:
#					print "UNCONNECTED"
#				return False
		
		# 1-lookahead rules:
		#
		# Border area size check
		#
		lhsneighs = len([ne for ne in n.GetAtomNeighbors() if ne in state.lhsborder])
		rhsneighs = len([ne for ne in m.GetAtomNeighbors() if ne in state.rhsborder])
		
		if (self.ISOMODE == GraphVF.SUBISOMORPH and lhsneighs < rhsneighs) or (self.ISOMODE == GraphVF.ISOMORPH and lhsneighs != rhsneighs):
			if self.verbose:
				print "BORDER SIZE MISMATCH (%d/%d)" % (lhsneighs, rhsneighs)
			return False
		
		# 2-lookahead rules:
		#
		# Remote area size check
		#
		num1 = len([ne for ne in n.GetAtomNeighbors() if ne not in state.lhsborder and ne not in state.map.keys()])
		num2 = len([ne for ne in m.GetAtomNeighbors() if ne not in state.rhsborder and ne not in state.map.values()])
		
		if (self.ISOMODE == GraphVF.SUBISOMORPH and num1 < num2) or (self.ISOMODE == GraphVF.ISOMORPH and num1 != num2):
			if self.verbose:
				print "UNMAPPED REGION SIZE MISMATCH (%d/%d)" % (num1, num2)
			return False
		
		return True
	
	
	def match(self, state):
		if self.boolmode:
			if self.isofound:
				return
		
		
		self.itercount += 1
		
		if len(state) == len(self.graph2.atoms):
			if len(self.graph1.atoms) == len(self.graph2.atoms):
				if self.verbose:
					print "isomorphism found:", state
			else:
				if self.verbose:
					print "subgraph isomorphism found:", state
			self.solutions.append(state)
			self.isofound = True
		else:
			if self.verbose:
				print " "*len(state), "state", state #, "[", lhs, rhs, "]"
			
			for n,m in self.candidatepairs(state):
				if self.verbose:
					print " "*len(state), "cp [%s] [%s]" % (n,m) #state
				if self.feasible(state, n, m):
					
					newstate = Statenode()
					newstate.map = state.map.copy()
					newstate.map[n] = m
					
					newstate.lhsborder = state.lhsborder[:]
					try:
						newstate.lhsborder.remove(n)
					except:
						pass
					newstate.rhsborder = state.rhsborder[:]
					try:
						newstate.rhsborder.remove(m)
					except:
						pass
					
					for ne in n.GetAtomNeighbors():
						if ne not in newstate.lhsborder and ne not in newstate.map.keys():
							newstate.lhsborder.append(ne)
					for ne in m.GetAtomNeighbors():
						if ne not in newstate.rhsborder and ne not in newstate.map.values():
							newstate.rhsborder.append(ne)
					
					newstate.lhsborder.sort()
					newstate.rhsborder.sort()
					
					self.match(newstate)
	
	def is_isomorphic(self):
		self.ISOMODE = ISOMORPH
		
		s = Statenode()
		s.lhsborder = []
		s.rhsborder = []
		s.map = {}
		
		self.boolmode = True
		self.match(s)
		
		return self.isofound
		
	
	def isomorph(self):
		self.ISOMODE = ISOMORPH
		
		s = Statenode()
		s.lhsborder = []
		s.rhsborder = []
		s.map = {}
		
		self.match(s)
		
		return self.solutions
	
	def subisomorph(self):
		self.ISOMODE = SUBISOMORPH
		
		s = Statenode()
		s.lhsborder = []
		s.rhsborder = []
		s.map = {}
		
		self.match(s)
		
		return self.solutions





# vf isomorphism for two sets of atoms/bonds, doesn't need graph or molecule structure,
# only atoms/bonds
class SubRegionVF:
	def __init__(self, atoms1, bonds1, atoms2, bonds2):
		self.atoms1 = atoms1
		self.atoms2 = atoms2
		self.bonds1 = bonds1
		self.bonds2 = bonds2
		self.itercount = 0
		self.solutions = []
		self.verbose = None
		
#		print "atoms1"
#		for a in self.atoms1:
#			print a, ",".join(map(str,a.neighs))
#		print "atoms2"
#		for a in self.atoms2:
#			print a, ",".join(map(str,a.neighs))
		
	
	# candidates that apply for next pair
	def candidatepairs(self, state):
		if state.lhsborder and state.rhsborder:
			rhs = state.rhsborder[0]
			for lhs in state.lhsborder:
				yield lhs,rhs
		
		# borders are empty!
		elif not state.lhsborder and not state.rhsborder:
			m = min( set(self.atoms2) - set(state.map.values()) )
			for n in self.atoms1:
				if n not in state.map.keys():
					yield n,m


	def feasible(self, state, n, m):
		#  atoms have to match
		if n.symbol != m.symbol:
			if self.verbose:
				print "ATOM MISMATCH"
			return False
		
		# connectedness-rule:
		# see that M[neighbors(n)] <=> neighbors(m)
		# i.e. that n's neighbors in M are m's neighbors
		#
		n_mapped_neighs = [ ne for ne in n.GetAtomNeighbors() if ne in state.map.keys()   and ne.GetBond(n) in self.bonds1]
		m_mapped_neighs = [ ne for ne in m.GetAtomNeighbors() if ne in state.map.values() and ne.GetBond(m) in self.bonds2]
		
		if len(n_mapped_neighs) != len(m_mapped_neighs):
			if self.verbose:
				print "MAPPED AREA SIZE MISMATCH"
			return False
		
		for node in n_mapped_neighs:
			if not state.map[node] in m_mapped_neighs:
				if self.verbose:
					print "MAPPED AREAS DIFFER"
				return False
		
		# border size rule:
		# size of border has to match
		#
		lhsneighs = len([ne for ne in n.GetAtomNeighbors() if ne in state.lhsborder and ne.GetBond(n) in self.bonds1])
		rhsneighs = len([ne for ne in m.GetAtomNeighbors() if ne in state.rhsborder and ne.GetBond(m) in self.bonds2])
		
		if lhsneighs != rhsneighs:
			if self.verbose:
				print "BORDER SIZE MISMATCH"
			return False
		
		# 2-lookahead rules:
		#
		# number of atoms have to match, which are in remote regions
		#  (MOLECULE \ [MAPPED_AREAS U BORDER]) 
		#
		num1 = 0
		num2 = 0
		for ne in n.GetAtomNeighbors():
			if ne in self.atoms1 and ne.GetBond(n) in self.bonds1 and ne not in state.lhsborder and ne not in state.map.keys():
				num1 += 1
		for ne in m.GetAtomNeighbors():
			if ne in self.atoms2 and ne.GetBond(m) in self.bonds2 and ne not in state.rhsborder and ne not in state.map.values():
				num2 += 1
		
		if num1 != num2:
			if self.verbose:
				print "REMOTE AREA SIZE MISMATCH"
			return False
		
		return True
	

	def match(self, state):
		self.itercount += 1
		
		if len(state) == len(self.atoms2):
			if self.verbose:
				print "isomorphism found:", state
			self.solutions.append(state.map)
		else:
			if self.verbose:
				print " "*len(state), "state", state#, "[", lhs, rhs, "]"
			
			for n,m in self.candidatepairs(state):
				if self.verbose:
					print " "*len(state), "cp", n, m, #state
				if self.feasible(state, n, m):
					
					if self.verbose:
						print
					
					newstate = Statenode()
					newstate.map = state.map.copy()
					newstate.map[n] = m
					
					newstate.lhsborder = state.lhsborder[:]
					try:
						newstate.lhsborder.remove(n)
					except:
						pass
					newstate.rhsborder = state.rhsborder[:]
					try:
						newstate.rhsborder.remove(m)
					except:
						pass
					
					for ne in n.GetAtomNeighbors():
						if ne in self.atoms1 and ne.GetBond(n) in self.bonds1 and ne not in newstate.lhsborder and ne not in newstate.map.keys():
							newstate.lhsborder.append(ne)
					for ne in m.GetAtomNeighbors():
						if ne in self.atoms2 and ne.GetBond(m) in self.bonds2 and ne not in newstate.rhsborder and ne not in newstate.map.values():
							newstate.rhsborder.append(ne)
					
					newstate.lhsborder.sort()
					newstate.rhsborder.sort()
					
					self.match(newstate)
	
	def isomorph(self):
		s = Statenode()
		s.lhsborder = []
		s.rhsborder = []
		s.map = {}
		
		self.match(s)
		
		return self.solutions





class MoleculeVF:
	ISOMORPH = 0
	SUBISOMORPH = 1
	MAPPINGISOMORPH = 2
	
	def __init__(self, mol1, mol2):
		self.mol1 = mol1
		self.mol2 = mol2
		self.itercount = 0
		self.ISOMODE = None
		self.solutions = []
		self.verbose = None
	
	def candidatepairs(self, state):
		if state.lhsborder and state.rhsborder:
			rhs = state.rhsborder[0]
			for lhs in state.lhsborder:
				yield lhs,rhs
		
		# borders are empty!
		elif not state.lhsborder and not state.rhsborder:
			m = min( set(self.mol2.graph.atoms) - set(state.map.values()) )
			for n in self.mol1.graph.atoms:
				if n not in state.map.keys():
					yield n,m


	def feasible(self, state, n, m):
		# the regions of G1 and G2 are divided into three areas:
		# - mapped regions
		# - adjacent border of the mapped region
		# - remote regions
		#
		# basically we have to check that the size of remote and border regions
		# has to match. also 
		#	
		
		# semantic rule:
		#  atoms have to match
		#
		if n.symbol != m.symbol:
			if self.verbose:
				print "ATOM MISMATCH"
			return False
		
		# connectedness-rule
		#
		# see that M[neighbors(n)] == neighbors(m)
		# i.e. that n's neighbors in M are m's neighbors
		#
		
		n_mapped_neighs = [ ne for ne in n.GetAtomNeighbors() if ne in state.map.keys() ]
		m_mapped_neighs = [ ne for ne in m.GetAtomNeighbors() if ne in state.map.values() ]
		
		if len(n_mapped_neighs) != len(m_mapped_neighs):
			if self.verbose:
				print "UNCONNECTED"
			return False
		
		for node in n_mapped_neighs:
			if not state.map[node] in m_mapped_neighs:
				if self.verbose:
					print "UNCONNECTED"
				return False
		
		# 0-lookahead rules:
		#  in the partial mapping:
		#   predecessors of candidate pair have to be paired together
		#   successors of candidate pair have to be paired together
		#
		# i.e. if candidate pair is (3,7), in the partial mapping
		#  it is not allowed to have a pair (5,3)
		#
	#	for nn,mm in state.map.items():
	#		if not ((nn < n and mm < m) or (nn > n and mm > m)):
	#			if options.verbose:
	#				print "RULE 0"
	#			return False
		
		
		# 1-lookahead rules:
		#
		# the number of predecessors of candidate pair belonging to the "border" 
		#  have to match
		# the number of predecessors of candidate pair belonging to the "border" 
		#  have to match
		#
		# this means, that if candidate pair is (3,7), the number of nodes [1..2] and
		#  [1...6] which lie at the 'border' of the mapping has to be same
		#
		# border means the immediate nodes which are not in the mapping, thus
		#  only unmapped regions have to processed
		#
		
		
		lhsneighs = 0
		rhsneighs = 0
		for ne in n.GetAtomNeighbors():
			if ne in state.lhsborder:
				lhsneighs += 1
		for ne in m.GetAtomNeighbors():
			if ne in state.rhsborder:
				rhsneighs += 1
		
		if (self.ISOMODE == MoleculeVF.SUBISOMORPH and lhsneighs < rhsneighs) or (self.ISOMODE == MoleculeVF.ISOMORPH and lhsneighs != rhsneighs):
			return False
		
		# 2-lookahead rules:
		#
		# number of predecessors have to match, which are in
		#  (MOLECULE - MAPPED_AREAS - BORDER) 
		#
		# basically the ares which are not immediately connected have to have 
		#  same number of pred and succ nodes
		#
		
		
		num1 = 0
		num2 = 0
		for ne in n.GetAtomNeighbors():
			if ne not in state.lhsborder and ne not in state.map.keys():
				num1 += 1
		for ne in m.GetAtomNeighbors():
			if ne not in state.rhsborder and ne not in state.map.values():
				num2 += 1
		
		if (self.ISOMODE == MoleculeVF.SUBISOMORPH and num1 < num2) or (self.ISOMODE == MoleculeVF.ISOMORPH and num1 != num2):
			return False
		
		return True
	

	def match(self, state):
		self.itercount += 1
		
		if len(state) == len(self.mol2.graph.atoms):
			if len(self.mol1.graph.atoms) == len(self.mol2.graph.atoms):
				if self.verbose:
					print "isomorphism found:", state
			else:
				if self.verbose:
					print "subgraph isomorphism found:", state
			self.solutions.append(state)
		else:
			if self.verbose:
				print " "*len(state), "state", state#, "[", lhs, rhs, "]"
			
			for n,m in self.candidatepairs(state):
				if self.verbose:
					print " "*len(state), "cp", cp[0], cp[1] #state
				if self.feasible(state, n, m):
					
					newstate = Statenode()
					newstate.map = state.map.copy()
					newstate.map[n] = m
					
					newstate.lhsborder = state.lhsborder[:]
					try:
						newstate.lhsborder.remove(n)
					except:
						pass
					newstate.rhsborder = state.rhsborder[:]
					try:
						newstate.rhsborder.remove(m)
					except:
						pass
					
					for ne in n.GetAtomNeighbors():
						if ne not in newstate.lhsborder and ne not in newstate.map.keys():
							newstate.lhsborder.append(ne)
					for ne in m.GetAtomNeighbors():
						if ne not in newstate.rhsborder and ne not in newstate.map.values():
							newstate.rhsborder.append(ne)
					
					newstate.lhsborder.sort()
					newstate.rhsborder.sort()
					
					self.match(newstate)


	def isomorph(self):
		self.ISOMODE = MoleculeVF.ISOMORPH
		
		s = Statenode()
		s.lhsborder = []
		s.rhsborder = []
		s.map = {}
		
		self.match(s)
		
		return self.solutions
	
	def subisomorph(self):
		self.ISOMODE = MoleculeVF.SUBISOMORPH
		
		s = Statenode()
		s.lhsborder = []
		s.rhsborder = []
		s.map = {}
		
		self.match(s)
		
		return self.solutions


class MapPair:
	def __init__(self, left, right, mapping = None):
		self.lhs = left
		self.rhs = right
		self.mapping = mapping
	
	def symbol(self):
		return self.lhs.symbol
	
	def neighbors(self):
		neighs = set()
#		print "finding neighbors of", self
		for ne in self.lhs.GetAtomNeighbors():
			neighs.add( MapPair(ne, self.mapping.mapping[ne], self.mapping) )
#			print ne, MapPair(ne, self.mapping.mapping[ne], self.mapping) 
		for ne in self.rhs.GetAtomNeighbors():
#			if ne not in neighs:
			for a in self.mapping.mapping.keys():
				if self.mapping.mapping[a] is ne:
#					print ne
					neighs.add( MapPair(a, ne, self.mapping) )
#					print ne, MapPair(a, ne, self.mapping) 
		
		return neighs
	
	# ==
	def __eq__(self, other):
		return self.lhs == other.lhs and self.rhs == other.rhs
	
	# !=
	def __ne__(self,other):
		return not self == other
	
	# TODO hash-funktio ei oikeen futaa
	def __hash__(self):
#		print hash(self.lhs) + hash(self.rhs)
		return hash(self.lhs) * hash(self.rhs)
	
	def __str__(self):
#		if self.lhs.graph.molecule.mol_id == "C03961":
		return "%d-%d" % (self.lhs.id, self.rhs.id)
#		else:
#			return "%d-%d" % (self.rhs.id, self.lhs.id)
#		return "%s,%s" % (self.lhs, self.rhs)
	
	
	# self > other
	def __gt__(self, other):
		if self.lhs.id > other.lhs.id:
			return True
		elif self.lhs.id < other.lhs.id:
			return False
		else:
			if self.rhs.id > other.rhs.id:
				return True
			else:
				return False
		
		
	# self < other
	def __lt__(self, other):
		if self.lhs.id < other.lhs.id:
			return True
		elif self.lhs.id > other.lhs.id:
			return False
		else:
			if self.rhs.id < other.rhs.id:
				return True
			else:
				return False



class MappingVF:
	def __init__(self, map1, map2):
		self.itercount = 0
		self.solutions = []
		self.verbose = None
		
		self.map1 = map1
		self.map2 = map2
		
		self.map1pairs = set()
		for pair in self.map1.mapping.items():
#			print pair[0], pair[1]
			if pair[0].graph.molecule in pair[0].graph.molecule.reaction.subs:
				self.map1pairs.add( MapPair(pair[0],pair[1],self.map1) )
			else:
				self.map1pairs.add( MapPair(pair[1],pair[0],self.map1) )
#		print
		self.map2pairs = set()
		for pair in self.map2.mapping.items():
#			print pair[0], pair[1]
#			if MapPair(pair[0],pair[1],self.map2) in self.map2pairs:
#				print "jepujee, on jo"
			if pair[0].graph.molecule in pair[0].graph.molecule.reaction.subs:
				self.map2pairs.add( MapPair(pair[0],pair[1],self.map2) )
			else:
				self.map2pairs.add( MapPair(pair[1],pair[0],self.map2) )
		
	
	def candidatepairs(self, state):
		if state.lhsborder and state.rhsborder:
			rhs = min(state.rhsborder)
			for lhs in state.lhsborder:
				yield lhs,rhs
		
		# borders are empty!
		elif not state.lhsborder and not state.rhsborder:
			m = min( self.map2pairs - set(state.map.values()) )
			for n in self.map1pairs - set(state.map.keys()):
#				if n not in state.map.keys():
				yield n,m
	
	
	# 'n' and 'm' are MapPair's
	def feasible(self, state, n, m):
#		print "testing feasibility for", n, m
		# ATOM COMPATIBILITY
		if n.symbol() != m.symbol():
#			print "ATOM FAIL", n.symbol(), m.symbol()
			return False
		
		# CONNECTEDNESS
		# M(mapped_neighbors(n)) == mapped_neighbors(m)
		n_neighbors_mapped = set( map( state.map.__getitem__, n.neighbors() & set(state.map.keys()) ) )
		m_neighbors = m.neighbors() & set(state.map.values())
		
		if n_neighbors_mapped != m_neighbors:
#			print "CONNECTEDNESS FAIL", n_neighbors_mapped, m_neighbors
			return False
		
		# BORDER SIZE MATCHING
		n_neighs_border = state.lhsborder & n.neighbors()
		m_neighs_border = state.rhsborder & m.neighbors()
		
		if len(n_neighs_border) != len(m_neighs_border):
#			print "BORDER FAIL", n_neighs_border, m_neighs_border
			return False
		
		# REMOTE SIZE MATCHING
		n_remote = self.map1pairs - set(state.map.keys())   - state.lhsborder
		m_remote = self.map2pairs - set(state.map.values()) - state.rhsborder
		
#		print state.lhsborder, state.rhsborder
#		print len(self.map1pairs)
#		print len(self.map2pairs)
#		print state.map
		
		if len(n_remote) != len(m_remote):
#			print "REMOTE FAIL", len(n_remote), len(m_remote)
			return False
		
		return True
	
	def match(self, state):
		self.itercount += 1
#		if self.itercount > 5:
#			sys.exit(1)
		if len(state) == len(self.map2pairs):
			if self.verbose:
				print "isomorphism found:"#, state
			self.solutions.append(state)
		else:
			if self.verbose:
				print " "*len(state), "state", state#, "[", lhs, rhs, "]"
			
			for n,m in self.candidatepairs(state):
				if self.verbose:
					print " "*(len(state)+1), "cp", n, m #state
				if self.feasible(state, n, m):
					
					newstate = Statenode()
					newstate.map = state.map.copy()
					newstate.map[n] = m
					
#					print "building new borders"
					
					newstate.lhsborder = state.lhsborder.copy()
					newstate.lhsborder.discard(n)
					newstate.rhsborder = state.rhsborder.copy()
					newstate.rhsborder.discard(m)
					
					for ne in n.neighbors():
						if ne not in newstate.lhsborder and ne not in newstate.map.keys():
							newstate.lhsborder.add(ne)
					for ne in m.neighbors():
						if ne not in newstate.rhsborder and ne not in newstate.map.values():
							newstate.rhsborder.add(ne)
					
					self.match(newstate)
	
	
	def isomorph(self):
		s = Statenode()
		
		self.match(s)
		
		return self.solutions
	




# tarviin algon joka ottaa kaksi reaktioverkkoa, ja laskee (sub)isomorfian.
# reaktioverkon kaaret annotoidaan: 0 ei muutosta, -1 katkos, +1 uusi kaari
# solmut ovat atomit
# lis�ksi reaktioverkot ovat superverkkoja: kullakin solmulla on lista 
# reaktioista tai compoundeista jotka kuuluvat siihen (tai jopa C00003:7)
# kaarilla on reaktiokohtainen katkostieto (k�tev�sti dict)
#


#class ReactionGraphVF(VF):
#	
#	
#	comes from VF
#	def candidatepairs(self, state)
#	
#	comes from VF
#	def feasible(self, state, n, m)
	






class Statenode:
	def __init__(self):
		self.map = {}
		self.lhsborder = set()
		self.rhsborder = set()
	
	def __len__(self):
		return len(self.map)
	
	def __str__(self):
		s = "{"
		for n,m in self.map.items():
			s += " %s|%s" % (n.id, m.id)
		s = s.strip("\n")
		s += " } "
		
		return s + "[" + ",".join(map(str, self.lhsborder)) +"] ["+ ",".join(map(str,self.rhsborder)) + "]"






if __name__ == "__main__":
	global options
	
	parser = optparse.OptionParser()
	parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose mode")
	parser.add_option("-q", "--quiet", dest="quiet", action="store_true", default=False, help="Quiet mode")
	
	(options, args) = parser.parse_args()
	
	if len(args) not in [1,2]:
		parser.error("give at least one molecule")
	
	if len(args) == 1:
		args.append(args[0])
	
	
	mol1 = Kegg.Molecule(args[0])
	mol2 = Kegg.Molecule(args[1])
	
	if len(mol1.graph.atoms) < len(mol2.graph.atoms):
		mol1,mol2 = mol2,mol1
	
	vf = VF(mol1,mol2)
	vf.verbose = options.verbose
	
	if len(mol1.graph.atoms) > len(mol2.graph.atoms):
		print "starting vf subgraph isomorphism algorithm for molecules %s (%d atoms) and %s (%d atoms) ..." % (args[0],len(mol1.graph.atoms), args[1], len(mol2.graph.atoms))
		
		starttime = time.time()
		sols = vf.subisomorph()
		endtime = time.time()
	else:
		print "starting vf isomorphism algorithm for molecules %s and %s of size %d atoms ..." % (args[0], args[1], len(mol1.graph.atoms))

		starttime = time.time()
		sols = vf.isomorph()
		endtime = time.time()
	
	
	print "Solutions:", len(sols)
	print "Iterations:", vf.itercount
	print "%.2f secs" % (endtime - starttime)



