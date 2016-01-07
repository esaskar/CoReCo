# implements undirected graph 
#
# basic philosophy is that you have to create the nodes and edges yourself, then add them to graph
# by AddNode/Edge methods
#



class Graph:
	idcounter = 0
	
	def __init__(self):
		self.id = Graph.idcounter
		Graph.idcounter += 1
		
		Node.idcounter = 0 # reset node and edge ids
		Edge.idcounter = 0
		
		self.nodes = []
		self.edges = []
	
	# adds a node
	def addnode(self):
		n = Node(self)
		self.nodes.append(n)
		return n
	
	# adds a node to existing node
	def addnodeto(self, to):
		if not to:
			return
		
		n = Node(self)
		e = Edge(n,to, self)
		to.addedge(e)
		n.addedge(e)
		
		self.nodes.append(n)
		self.edges.append(e)
		return n
	
	def removenode(self, n):
		while n.edges:
			self.removeedge(n.edges[0])
		self.nodes.remove(n)
		
	
#	def addedge(self):
#		e = Edge(self)
#		self.edges.append(e)
#		return e
	
	def addedgeto(self, n1, n2):
		if not n1 or not n2:
			return
		
		if n1.getedge(n2):  # already exists
			return
		
		e = Edge(n1, n2, self)
		
		e.source.addedge(e)
		e.target.addedge(e)
		
		self.edges.append(e)
		return e
	
	def removeedge(self, p1, p2=None):
		if p2:
			e = p1.getedge(p2)
		else:
			e = p1
		
		if not e:  # doesn't exist
			return
		
#		if e.source.id == 1 or e.target.id == 1:
#			print map(str, e.source.edges), map(str, e.source.neighs)
		
		e.source.edges.remove(e)
		e.target.edges.remove(e)
		e.source.neighs.remove(e.target)
		e.target.neighs.remove(e.source)
		self.edges.remove(e)
	
	# for a set of nodes, get their border
	def neighbors(self, nodes):
		return list(set([ne for n in nodes for ne in n.neighbors() if ne not in nodes]))
		
	
	
	def __str__(self):
		s =  "graph %d (%d/%d)\n" % (self.id, len(self.nodes), len(self.edges))
		for n in self.nodes:
			s += str(n) + "\n"
		for e in self.edges:
			s += str(e) + "\n"
		return s.strip()


class Node:
	idcounter = 0
	
	def __init__(self, graph = None):
		self.id = Node.idcounter
		Node.idcounter += 1
		
		self.graph = graph
		self.edges = []
		self.neighs = []
	
	
	def addedge(self, e):
		if not e:
			return
		
		self.edges.append(e)
		self.neighs.append( e.source if e.target is self else e.target )
	
	def degree(self):
		return len(self.edges)
	
	def neighbors(self):
		return self.neighs
	
	def connected(self, other):
		if other in self.neighs:
			return True
		return False
	
	def edges(self):
		return self.edges
	
	def getedge(self, node):
		for e in self.edges:
			if e.source is node or e.target is node:
				return e
		return None
	
	def getnode(self, edge):
		return edge.source if edge.target is self else edge.target
	
	def __str__(self):
		return "node %d" % (self.id)


class Edge:
	idcounter = 0
	
	def __init__(self, n1, n2, graph = None):
		self.id = Edge.idcounter
		Edge.idcounter += 1
		
		self.graph = graph
		self.source = n1
		self.target = n2
	
	def other(self, n):
		return self.source if n == self.target else self.target
	
	def __str__(self):
		return "edge %d (%d <-> %d)" % (self.id, self.source.id, self.target.id)
	







