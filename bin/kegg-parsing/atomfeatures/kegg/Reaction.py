# -*- coding: iso-8859-15 -*-

from Molecule import *
from Mapping import *
from ReactionGraph import *
import Kegg

class Reaction:
	def __init__(self, reaction_id = None, block = None):
		self.reaction_id = reaction_id
		self.ligand = self.reaction_id
		self.equation = None #str
		self.subs = []
		self.prods = []
		self.reactiongraph = None
		
		self.height = 0
		self.width = 0
		
		self.mappings = []
#		self.arita_mapping = None
#		self.kegg_mapping = None
#		self.astar_mapping = None
		self.rpairs = []
		
		if reaction_id:
			self.parse_reaction(block)
			self.read_reactants()
	#		self.compute_mapping_colors()
		
	
	def parse_reaction(self, block = None):
		
		data = block
		
		if not block:
			f = open(Kegg.Kegg.REACTION_FILE)
			data = f.readlines()
			f.close()
		
		self.equation = ""
		lastkeyword = None
		reac = False
		for line in data:
			line = line.rstrip()
			
			if not line:
				continue
			
			if not line.startswith(" "):
				lastkeyword = line.split()[0]
			
			if line.startswith("ENTRY") and self.reaction_id in line:
				reac = True
			elif "///" in line and reac == True:
				break
			elif reac and ("EQUATION" in line or (line.startswith(" ") and lastkeyword == "EQUATION")):
				self.equation += line[12:].strip() + " "
				eq = True
			elif lastkeyword == "RPAIR" and reac:
				self.rpairs.append(line[12:].split()[1])
		
		words = self.equation.strip().split()
		self.subs = []
		self.prods = []
		
		ptr = self.subs
		
		for i in range(len(words)):
			w = words[i]
			
			if w == "+":
				continue
			elif w.isdigit():
				wnext = words[i+1]
				if re.search(r'C[0-9]{5}', wnext):
					for k in range(0,int(w)-1):
						wnext = re.findall(r'C[0-9]{5}', wnext)[0]
						ptr.append(wnext)
			elif w == "<=>":
				ptr = self.prods
			elif re.search(r'C[0-9]{5}', w):
				w = re.findall(r'C[0-9]{5}', w)[0]
				ptr.append(w)
		
	def read_reactants(self):
		left = 50
		max_right = 0
		top = 50
		max_bottom = 0
		
		for i in range(len(self.subs)):
			m = Molecule(self.subs[i], self, left, top)
			self.subs[i] = m
			top = m.max_y + 50
			if max_right < m.max_x:
				max_right = m.max_x
			if max_bottom < m.max_y:
				max_bottom = m.max_y
		
		left = max_right + 100
		top = 50
		
		for i in range(len(self.prods)):
			m = Molecule(self.prods[i], self, left, top)
			self.prods[i] = m
			top = m.max_y + 50
			if max_right < m.max_x:
				max_right = m.max_x
			if max_bottom < m.max_y:
				max_bottom = m.max_y
		
		self.subs = [m for m in self.subs if m.graph and m.graph.atoms] # filter atomless molecules
		self.prods = [m for m in self.prods if m.graph and m.graph.atoms] # filter atomless molecules
		
		# trimmed formula containing only non-empty molecules
		self.trim_formula = " + ".join([m.ligand for m in self.subs]) + " <=> " + " + ".join([m.ligand for m in self.prods])
		
		self.height = max_bottom + 50
		self.width = max_right + 50
		
	
	def read_mapping(self, t = "ASTAR", fname = None):
		if t == "ASTAR":
			self.mappings.append( Mapping(self) )
			self.mappings[-1].read_astar_mapping(fname)
		elif t == "KEGG":
			self.mappings.append( Mapping(self) )
			self.mappings[-1].read_kegg_mapping(fname)
		elif t == "ARITA":
			self.mappings.append( Mapping(self) )
			self.mappings[-1].read_arita_mapping(fname)
	
	def read_reactiongraph(self, mapfile = None):
		self.reactiongraph = ReactionGraph(self)
		self.reactiongraph.read(mapfile)
	
	
	def compute_mapping_colors(self):
		pass
	

	def write_header(self, f):
		f.write("<?xml version=\"1.0\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n\n<svg fill-opacity=\"1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" color-rendering=\"auto\" color-interpolation=\"auto\" stroke=\"black\" text-rendering=\"auto\" stroke-linecap=\"square\" width=\"%d\" stroke-miterlimit=\"10\" stroke-opacity=\"1\" shape-rendering=\"auto\" fill=\"black\" stroke-dasharray=\"none\" font-weight=\"normal\" stroke-width=\"1\" height=\"%d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"&apos;Dialog&apos;\" font-style=\"normal\" stroke-linejoin=\"miter\" font-size=\"12\" stroke-dashoffset=\"0\" image-rendering=\"auto\">\n" % (self.width, self.height))
		f.write("  <g>\n")
		
	def write_end(self, f):
		f.write("  </g>\n")
		f.write("</svg>\n")
	
	def write_background(self, f):
		f.write('    <g fill="rgb(238,238,238)" stroke="rgb(238,238,238)">\n      <rect width="%d" x="0" height="%d" y="0" stroke="none" />\n    </g>\n' % (self.width, self.height))
	
	def write_svg(self):
		f = open("./" + self.reaction_id + ".svg", "w")
		
		self.write_header(f)
		
		for m in self.subs + self.prods:
			m.write_svg(f)
			f.write("\n")
		
		self.write_end(f)
		
		f.close()
	
	def atoms(self):
		return [a for m in self.subs+self.prods for a in m.graph.atoms]
	
	def get_atoms(self, mol_id, atom_id):
		mol = None
		atom = None
		
		try:
			mols = [m for m in self.subs + self.prods if m.mol_id == mol_id]
		except:
			print "wanted mol", mol_id, ". have:", self.equation
		
		try:
			atoms = [a for m in mols for a in m.graph.atoms if a.id == int(atom_id)]
		except:
			print "wanted atom", atom_id, "have:"
			for a in mols[0].graph.atoms:
				print a.id
			return
		
		return atoms
	
	# checks whether atom numbering in rpair and mol-file are same
#	def atom_numbering_consistent(self):
#		for mol in self.subs + self.prods:
#			# read from rpair!
#			Mapping.read_rpair_data() 
	
	
	
	def has_hydrogens(self):
		for c in self.subs+self.prods:
			for a in c.graph.atoms:
				if a.symbol == "H" or a.symbol == "H+":
					return True
		
		return False
	
	def is_balanced(self):
		subatoms = sorted("".join([a.symbol for m in self.subs for a in m.graph.atoms]))
		prodatoms = sorted("".join([a.symbol for m in self.prods for a in m.graph.atoms]))
		
		return subatoms == prodatoms
	
	
	def __str__(self):
		return self.reaction_id
	
	
	
	
	
	
	
	
	
		