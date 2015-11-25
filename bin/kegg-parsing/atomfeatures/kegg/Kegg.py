import re, glob

from Reaction import *
from Molecule import *


class Kegg:
	
	MOL_FOLDER = "../../../data/Kegg/mol/"
	REACTION_FILE = "../../../data/Kegg/reaction"
	
	def __init__(self):
		self.reactions = []
		self.molecules = []
		
#		self.read_molecules()
		self.read_reactions()
	
	
	def read_molecules(self):
		files = glob.glob(Kegg.MOL_FOLDER + "*.mol")
		files = map(lambda fn: fn[fn.rfind("/"):-4], files) # take only short names
		
		for fn in files:
			m = Molecule(fn)
			self.molecules.append( m )
	
	
	def read_reactions(self):
		f = open(Kegg.REACTION_FILE)
		data = f.read()
		f.close()
		
		print "reaction file read into memory"
		
		blocks = data.split("\n///\n")
		blocks = [ b.split("\n") for b in blocks ]
		blocks = [ [x for x in b if x.strip() != ""] for b in blocks ]
		blocks = [ b for b in blocks if b ]
		
		print "reactions digested into blocks"
		
		blocks = blocks[0:100]
		
		for block in blocks:
#			print block[0]
			reac_id = None
			try:
				reac_id = re.findall(r'(R\d{5})', block[0])[0]
			except:
				print block, "failed"
			
			reac = Reaction(reac_id, block)
			
			if reac:
				self.reactions.append( reac )
			else:
				print "reaction failed"
				
			if len(self.reactions) % 100 == 0:
				print len(self.reactions), "reactions processed"
		
#		print "1000 reactions read"
	# read all the mappings
	# assumes reactions read
	def read_mappings(self):
		
#		rpairdata = Mapping.read_rpair()
#		
#		self.rpairs = rpairdata
#		
#		print "RPAIR data read into memory"
		
#		self.reactions[91].read_mapping(rpairdata)
		
		compdist = {}
		complete = 0
		
		for r in self.reactions:
			if r.is_defined():
				r.read_mapping("ASTAR")
				comp = r.mapping.completeness()
				if compdist.has_key(comp):
					compdist[comp] += 1
				else:
					compdist[comp] = 1
				if comp[0] == comp[1]:
					complete += 1
				print r, r.mappings[0]
#				print "%s: %d/%d" % (r.reaction_id, comp[0], comp[1])
#				print r, "mapped"
		
		print complete
		print "Mappings read"
		
		return compdist
	
	
	def __str__(self):
		return "%d reactions" % (len(self.reactions))













