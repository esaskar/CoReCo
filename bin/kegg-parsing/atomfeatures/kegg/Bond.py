# -*- coding: iso-8859-15 -*-

class Bond:
	bondcounter = 0
	
	def __init__(self, source, target, type = 1, id = bondcounter, graph = None):
		self.graph = graph
		self.source = source
		self.target = target
		self.type = int(type)
		self.id = id
	
	def other(self, atom):
		if self.source is atom:
			return self.target
		elif self.target is atom:
			return self.source
		return None
	
#	def __eq__(self, other):
#		return self.source is other.source and self.target is other.target and self.graph is other.graph
#	
#	def __ne__(self, other):
#		return not (self == other)
	
	
	def __str__(self):
		return "Bond " + str(self.id) + ": " + str(self.source) + " <-> " + str(self.target)
