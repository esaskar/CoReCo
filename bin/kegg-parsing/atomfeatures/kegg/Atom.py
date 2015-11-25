# -*- coding: iso-8859-15 -*-

from Bond import *
#from ReactionGraph import *


class Atom:
	atomcounter = 0

	def __init__(self, symbol, id = None, graph = None, x = None, y = None, color = None):
		self.graph = graph
		self.symbol = symbol
		self.bonds = []
		self.neighs = []
		
		self.x = x
		self.y = y
		self.color = color
		
		if id is None:
			self.id = Atom.atomcounter
		else:
			self.id = id
		
		Atom.atomcounter += 1
	
	def AddBond(self, bond):
		self.bonds.append(bond)
		if bond.source is self:
			self.neighs.append(bond.target)
		else:
			self.neighs.append(bond.source)

	def GetBond(self, atom):
		if not atom or atom is self:
			return None
		
		for b in self.bonds:
			if b.target is atom or b.source is atom:
				return b
		
		return None
	
	def GetAtomNeighbors(self):
		neighs = []
		
		for b in self.bonds:
			if b.source is not self:
				neighs.append(b.source)
			else:
				neighs.append(b.target)
		
		return neighs
	
	def GetBondNeighbors(self):
		return self.bonds
		
	def IsNeighbor(self, other):
		if other in self.neighs:
			return True
		return False
	
	def Degree(self):
		return len(self.bonds)
	
	def __str__(self):
		try:
			return "Atom %s (%d/%s)" % (self.symbol, self.id, self.graph.molecule.mol_id)
		except:
			return "Atom %s (%d)" % (self.symbol, self.id)
#		return str(self.id) + "/" + str(self.graph.molecule.mol_id)
	
	def __hash__(self):
		return id(self)
	
#	# self > other
#	def __gt__(self, other):
#		return self.id > other.id
#	
#	# self < other
#	def __lt__(self, other):
#		return self.id < other.id
#	
#	def __cmp__(self, other):
#		if self.id < other.id:
#			return -1
#		elif self.id > other.id:
#			return 1
#		return 0

#class SuperAtom(Atom):
#	def __init__(self, *arg, **kw):
#		Atom.__init__(self, arg, kw)
#		
#	def __str__(self):
#		return "Atom %s (%d/%s)" % (self.symbol, self.id, self.graph.reaction.reaction_id)




