from Reaction import *
from Molecule import *
from Graph import Graph

# metabolic network, a bi-partite
# each node has a type (reaction/molecule)
# each edge is a simple connector
# molecule nodes are not modeled yet...
class MetabolicNetwork():
	MOL_FOLDER = "/group/home/icomic/data/kegg/ligand/LATEST/mol/"
	REACTION_FILE = "/group/home/icomic/data/kegg/ligand/LATEST/reaction"
	
#	COFACTORS = set(["NAD+", "NADH", "ATP", "ADP", "AMP", "H2O", "GDP", "Diphosphate", "H+"])
	COFACTORS = set(["C00001", "C00002", "C00080", "C00008", "C00009", "C00007", "C00014", "C00020"])
	
	def __init__(self, reactionlist = None):
#		self.graph = Graph()
		self.reactions = []
		self.adjmatrix = []
		self.reactionlist = reactionlist
		
		self.read_reactions()
		self.create_network()
		
	
#	def read_molecules(self):
#		files = glob.glob(MetabolicNetwork.MOL_FOLDER + "*.mol")
#		files = map(lambda fn: fn[fn.rfind("/"):-4], files) # take only short names
#		
#		for fn in files:
#			m = Molecule(fn)
#			self.compounds.append( m )
	
	
	def rgneighs(self, reaction_id):
#		print "neighs of", reaction_id
#		for i in range(len(self.reactions)):
#			print self.reactions[i], self.adjmatrix[i]
		
		result = []
		
#		print reaction_id
		for i in range(len(self.reactions)):
#			print i, repr(self.reactions[i].ligand), repr(reaction_id), self.reactions[i].ligand == reaction_id
#			print i, "z"+ self.reactions[i].reaction_id + "z", reaction_id
			if self.reactions[i].reaction_id == reaction_id:
#				print "whazoo"
#				print i
				result = [self.reactions[j].reactiongraph for j in range(len(self.adjmatrix[i])) if self.adjmatrix[i][j] == 1]
				break
		
#		for r in result:
#			print r.reaction.ligand
			
#		print result
		return result
	
	
	def write_dot(self):
		# format is 
		# C00001
		# ...
		# R00001
		# ...
		# C00001 -> R00001
		# R00001 -> C00002
		
		f = open("megagraph.dot", "w")
		f.write("digraph g {\n")
		
		reactants = []
		for r in self.reactions:
			for reactant in r.subs + r.prods:
				reactants.append(reactant.mol_id)
		
		reactants = set(reactants)
		
		for r in reactants:
			f.write(r + ";\n")
		
		for r in self.reactions:
			f.write(r.reaction_id + ";\n")
		
		for r in self.reactions:
			for sub in r.subs:
				mol_id = sub.mol_id
				f.write(mol_id + " -> " + r.reaction_id + ";\n")
			for prod in r.prods:
				prod_id = prod.mol_id
				f.write(r.reaction_id + " -> " + prod_id + ";\n")
		
		f.write("}\n")
		f.close()
		
	
	
	def create_network(self):
		for r in self.reactions:
#			self.graph.AddNode(r)
#			r.reactants = set([m for m in r.subs + r.prods]
			r.reactant_ids = set([m.code for m in r.subs + r.prods])
		
		reacs = len(self.reactions)
		self.adjmatrix = [ [0 for i in range(reacs)] for x in range(reacs) ]
		
#		print "%d nodes created" % (len(self.graph.nodes))
		
		for i in range(reacs):
			count = 0
			
			for j in range(reacs):
				if i == j:
					continue
				
				r1 = self.reactions[i]
				r2 = self.reactions[j]
				
				if r1.reactant_ids.intersection(r2.reactant_ids):
#					self.graph.AddEdgeTo( r1,r2, ) #tuple(r1.reactants.intersection(r2.reactants)) )
					self.adjmatrix[i][j] = 1
					count += 1
#					print " edge between", str(r1), str(r2)
		
			print "node %d connected to %d nodes" % (i, count)
		
		for i in range(len(self.adjmatrix)-1,-1,-1):
			if sum(self.adjmatrix[i]) == 0:
				del self.adjmatrix[i]
				self.adjmatrix = [ row[:i]+row[i+1:] for row in self.adjmatrix]
				
		
	
	def read_reactions(self):
		f = open(MetabolicNetwork.REACTION_FILE)
		data = f.read()
		f.close()
		
		print "reaction file read into memory"
		
		blocks = data.split("\n///\n")
		blocks = [ b.split("\n") for b in blocks ]
		blocks = [ [x for x in b if x.strip() != ""] for b in blocks ]
		blocks = [ b for b in blocks if b ]
		
		print "reactions digested into blocks"
		
#		blocks = blocks[0:300]
		
		for block in blocks:
#			print block[0]
			reac_id = None
			try:
				reac_id = re.findall(r'(R\d{5})', block[0])[0]
			except:
				pass
#				print block, "failed"
			
			if not self.reactionlist or reac_id in self.reactionlist:
				reac = Reaction(reac_id, block)
				
				if reac:
					self.reactions.append( reac )
		#			else:
		#				print "reaction failed"
					
				if len(self.reactions) % 100 == 0:
					print len(self.reactions), "reactions processed"
			
	
	
	
	
	
	def size(self):
		return len(self.reactions)
	
	def __str__(self):
		return "Metabolic network with %d reactions" % (len(self.reactions))
	
	

if __name__ == "__main__":
	n = MetabolicNetwork()






