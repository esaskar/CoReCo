# -*- coding: iso-8859-15 -*-
#
# basic class representing molecules
# contains a 'graph' object which has the actual molecular structure
# if molecule doesn't have structure (e.g. "oxidizing agent" with no struc)
# then graph=None
# 
# even if no molecule couldn't be found, this object still is valid
#
# removes hydrogens by default
#
##

from MoleculeGraph import *
import Kegg
import math

class Molecule:
	def __init__(self, mol_id, reaction = None, left = None, top = None, molfile = None):
		self.mol_id = mol_id.strip()
		self.ligand = self.mol_id
#		self.code = mol_id[1:])
		self.code = self.mol_id
		self.graph = MoleculeGraph(self)
		self.reaction = reaction
		
		self.max_x = 0
		self.max_y = 0
		self.min_x = 0
		self.min_y = 0
		self.width = None
		self.height = None
		
		if molfile:
			self.read_molfile(molfile)
		else:
			self.read_mol()
		
		if self.graph:
			self.calibrate(left, top)
	
	def read_mol(self):
		self.read_molfile(Kegg.Kegg.MOL_FOLDER + self.mol_id + ".mol")
	
	def read_molfile(self, molfile):
		format = ""
		
#		print "opening file", self.ligand, molfile
		
		try:
			f = open(molfile)
		except IOError:
			self.graph = None # couldn't open file -> graph=null
#			print "file", molfile, "not found"
			return
		
		data = f.read()
		f.close()
		
		if "TRIPOS" in data:
			self.read_mol2file(molfile)
		else:
			self.read_mol1file(molfile)
	
	
	def read_mol1file(self, molfile):
		index = 1    # index of atoms
		atoms = {}   # index -> atom_obj
		bondindex = 1
	
		try:
			f = open(molfile)
		except IOError:
			self.graph = None # again, no file -> graph=null
			print "Couldn't construct Molecule using molfile: %s" % (molfile)
			return
		
		# turhat header-rivit
		section = "HEADER"
		f.readline()
		f.readline()
		f.readline()
		f.readline()
		
		for line in f:
			if line.startswith("M") or line.startswith("A") or line.startswith("R") :
				## added "A" and "R" to handle a mol file with R#'s (biocyc 0301121217)
				##A   65
				##R2
				##A   64
				##R1
				section = "OTHER"
			elif len(line.strip().split()) in [6,7]:
				section = "BONDS"
			elif len(line.strip().split()) in [11,13,16]: ## Added 11 as additional option to allow parsing of biocyc .mol-files (Merja Oja, 29.10.2013)
				section = "ATOMS"
			
			if section == "ATOMS":
				words = line.strip().split()
				
				symbol = words[3]
				x = float(words[0])
				y = float(words[1])
				
				a = Atom(symbol, index, graph = self.graph, x = x, y = y, color = "black")
				self.graph.AddAtom(a)
				atoms[index] = a
				index += 1
				
			elif section == "BONDS":
				words = line.strip().split()
				source = atoms[int(line[0:3].strip())]
				target = atoms[int(line[3:6].strip())]
				type = line[6:9].strip()
				
				b = Bond(source, target, type, bondindex)
				bondindex += 1
				
				len(line.rstrip())/3
				
				b.rest = []
				for i in range(3, len(line.rstrip())/3):
					b.rest.append(line[i*3:i*3+3].strip())
				
				self.graph.AddBond(b)
		
		f.close()
		
		self.remove_hydrogens()
	
	
	def read_mol2file(self, molfile):
		index = 1    # index of atoms
		atoms = {}   # index -> atom_obj
		bondindex = 1
		
		try:
			f = open(molfile)
		except IOError:
			self.graph = None
			return
		
		lines = f.readlines()
		f.close()
		
		section = ""
		
		for line in lines:
			line = line.strip()
			
			if "<TRIPOS>MOLECULE" in line:
				section = "MOLECULE"
				continue
			elif "<TRIPOS>ATOM" in line:
				section = "ATOM"
				continue
			elif "<TRIPOS>BOND" in line:
				section = "BOND"
				continue
			
			if section == "MOLECULE":
				pass
			elif section == "ATOM":
				words = line.split()
				
#				id = words[0]
				symbolinfo = words[5]
#				symbol = words[5].split(".")[0]
				charge = words[8]
				x = float(words[2])
				y = float(words[3])
				
				a = Atom(symbolinfo, index, graph = self.graph, x=x, y=y)
#				a.symbolinfo = symbolinfo
				self.graph.AddAtom(a)
				atoms[index] = a
				index += 1
				
			
			elif section == "BOND":
				id, source, target, type = line.split()
				source = atoms[source]
				target = atoms[target]
				
				b = Bond(id, source, target, type)
				bondindex += 1
				
				self.graph.AddBond(b)
		
		self.remove_hydrogens()
	
	
	def remove_hydrogens(self):
		for i in range(len(self.graph.atoms)-1, -1, -1):
			if self.graph.atoms[i].symbol == "H" or self.graph.atoms[i].symbol == "H+":
				a = self.graph.atoms[i]
				if a.neighs:
					ne = a.neighs[0] # only one neighbor (if any)
					ne.neighs.remove(a)
					self.graph.bonds.remove(a.bonds[0])
					ne.bonds.remove(a.bonds[0])
				
				del self.graph.atoms[i]
		
		# compress id's to go from 1...n in atoms and 1...m in bonds
		self.compress_ids()
	
	# compress id's to go from 1...n in atoms and 1...m in bonds
	# does nothing if already satisfied
	def compress_ids(self):
		i = 1
		for a in self.graph.atoms:
			a.id = i
			i += 1
		i = 1
		for b in self.graph.bonds:
			b.id = i
			i += 1
	
	
	def strip_atominfo(self):
		for a in self.graph.atoms:
			a.symbol = a.symbol.split(".")[0]
	
	def full_bondinfo(self):
		# complete bond types from "1" to "C1N" style
		for b in self.graph.bonds:
			srcsymbol = b.source.symbol.split(".")[0]
			tgtsymbol = b.target.symbol.split(".")[0]
			
			if srcsymbol > tgtsymbol:
				srcsymbol, tgtsymbol = tgtsymbol, srcsymbol
			
			b.type = srcsymbol + str(b.type) + tgtsymbol
	
	
	def calibrate(self, left, top):
		###############################################################
		#
		# scale coordinates such that distance between closest atoms is exactly 28px
		#
		###############################################################
		
		if not left:
			left = 50
		if not top:
			top = 50
		
		# find min-distance between any two atoms
		min_dist = 1000000
		min_dist_i = None
		min_dist_j = None
		
		for i in range(0,len(self.graph.atoms)-1):
			for j in range(i+1, len(self.graph.atoms)):
				dist = math.sqrt((self.graph.atoms[i].x - self.graph.atoms[j].x)**2 + (self.graph.atoms[i].y-self.graph.atoms[j].y)**2)
				if 0.01 < dist < min_dist:
					min_dist = dist
					min_dist_a = self.graph.atoms[i]
					min_dist_b = self.graph.atoms[j]
		
		# now minimum distance between any two atoms is 'min_dist'
		# lets scale this value to, say, '10'
		# do this by multiplying all coordinates with '10/min_dist'
		
		for i in range(0,len(self.graph.atoms)):
			self.graph.atoms[i].x *= 32 / min_dist
			self.graph.atoms[i].y *= 32 / min_dist
		
		
		###################################################
		#
		# make all coordinates positive and left-top corner -most atom into (50,50) coordinate
		#
		###################################################
		
		max_x = -1000
		min_x = 1000
		max_y = -1000
		min_y = 1000
		
		for atom in self.graph.atoms:
			if atom.x < min_x:
				min_x = atom.x
			if atom.x > max_x:
				max_x = atom.x
			if atom.y < min_y:
				min_y = atom.y
			if atom.y > max_y:
				max_y = atom.y
		
		for i in range(len(self.graph.atoms)):
			atom = self.graph.atoms[i]
			self.graph.atoms[i].x = atom.x - min_x + left
			self.graph.atoms[i].y = atom.y - min_y + top
	
		for atom in self.graph.atoms:
			if atom.x < min_x:
				min_x = atom.x
			if atom.x > max_x:
				max_x = atom.x
			if atom.y < min_y:
				min_y = atom.y
			if atom.y > max_y:
				max_y = atom.y
	
		self.width = int(max_x - min_x)
		self.height = int(max_y - min_y)
			
		self.max_x = max_x
		self.max_y = max_y
		self.min_x = min_x
		self.min_y = min_y
	
	
	
	def write_header(self, f):
		f.write("    <g>\n")
	
	def write_end(self, f):
		f.write("    </g>\n")
	
	def write_background(self, f):
		f.write('      <g fill="rgb(238,238,238)" stroke="rgb(238,238,238)">\n      <rect width="%d" x="0" height="%d" y="0" stroke="none" />\n    </g>\n' % (self.width, self.height))
	
	
	def write_atoms(self, f):
		f.write("      <g font-size=\"15px\" fill=\"rgb(0,0,0)\" text-rendering=\"optimizeLegibility\" font-family=\"Helvetica\" shape-rendering=\"geometricPrecision\" stroke=\"rgb(0,0,0)\">\n")
		for atom in self.graph.atoms:
			f.write("        <rect x=\"%d\" y=\"%d\" fill=\"white\" width=\"15\" height=\"11\" stroke=\"none\" />\n" % (atom.x+5,atom.y-2))
			if atom[2] == "C":
				continue
			f.write("        <rect x=\"%d\" y=\"%d\" fill=\"white\" width=\"15\" height=\"16\" stroke=\"none\" />\n" % (atom.x-7,atom.y-8))
			
		f.write("      </g>\n")
		
		
		
		f.write("      <g font-size=\"15px\" fill=\"rgb(0,0,0)\" text-rendering=\"optimizeLegibility\" font-family=\"Helvetica\" shape-rendering=\"geometricPrecision\" stroke=\"rgb(0,0,0)\">\n")
		
		for atom in self.graph.atoms:
			f.write("        <text font-size=\"11px\" fill=\"blue\" xml:space=\"preserve\" x=\"%d\" y=\"%d\" stroke=\"none\">%s</text>\n" % (atom.x+6,atom.y+8,str( atom[3]+1 )))
			
			if atom[2] == "C":
				continue
			
			f.write("        <text fill=\"black\" xml:space=\"preserve\" x=\"%d\" y=\"%d\" stroke=\"none\">%s</text>\n" % (atom.x-5,atom.y+5,str(atom[2])))
		
		f.write("      </g>\n")
		
	def write_bonds(self, f):
		f.write("      <g text-rendering=\"optimizeLegibility\" stroke-width=\"1.1\" shape-rendering=\"geometricPrecision\">\n")
		
		for bond in self.graph.bonds:
#			print bond
			x1 = self.graph.atoms[bond[0]].x
			x2 = self.graph.atoms[bond[1]].x
			y1 = self.graph.atoms[bond[0]].y
			y2 = self.graph.atoms[bond[1]].y
			
			if bond[2] == 2:
				# get direction
				dir_x = (x2-x1) / 6
				dir_y = (y1-y2) / 6
				
				f.write("        <line y2=\"%f\" fill=\"none\" x1=\"%f\" x2=\"%f\" y1=\"%f\" />\n" % (y2+dir_x+dir_y,x1+dir_y+dir_x,x2+dir_y-dir_x,y1+dir_x-dir_y))
				f.write("        <line y2=\"%f\" fill=\"none\" x1=\"%f\" x2=\"%f\" y1=\"%f\" />\n" % (y2,x1,x2,y1))
			else:
				f.write("        <line y2=\"%f\" fill=\"none\" x1=\"%f\" x2=\"%f\" y1=\"%f\" />\n" % (y2,x1,x2,y1))
		
		f.write("      </g>\n")
	
	
	def write_svg(self, f):
		self.write_header(f)
#		self.write_background(f)
		self.write_bonds(f)
		self.write_atoms(f)
		self.write_end(f)
	
	def get_atom(self, id):
#		print "finding atom", id, "for mol", self.mol_id
		for a in self.graph.atoms:
			if a.id == id:
				return a
		
		return None
	
	def compute_borders(self):
		width = max([a.x for a in self.graph.atoms]) - min([a.x for a in self.graph.atoms])
		height = max([a.y for a in self.graph.atoms]) - min([a.y for a in self.graph.atoms])
		
		return width + 100, height + 100
	
	
	def svg(self, fname = None):
		if not fname:
			fname = self.ligand + ".svg"
		
		
		width, height = self.compute_borders()
		
		# header
		s = "<?xml version=\"1.0\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n\n<svg fill-opacity=\"1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" color-rendering=\"auto\" color-interpolation=\"auto\" stroke=\"black\" text-rendering=\"auto\" stroke-linecap=\"square\" width=\"%d\" stroke-miterlimit=\"10\" stroke-opacity=\"1\" shape-rendering=\"auto\" fill=\"black\" stroke-dasharray=\"none\" font-weight=\"normal\" stroke-width=\"1\" height=\"%d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"&apos;Dialog&apos;\" font-style=\"normal\" stroke-linejoin=\"miter\" font-size=\"12\" stroke-dashoffset=\"0\" image-rendering=\"auto\">\n  <g>\n" % (width, height)
		
		# bonds
		s += "      <g text-rendering=\"optimizeLegibility\" stroke-width=\"1.1\" shape-rendering=\"geometricPrecision\">\n"
		for b in self.graph.bonds:
			x1 = b.source.x
			x2 = b.target.x
			y1 = b.source.y
			y2 = b.target.y
			
			s += "        <line x1=\"%f\" fill=\"none\" y1=\"%f\" x2=\"%f\" y2=\"%f\" />\n" % (x1,y1,x2,y2)
		s += "      </g>\n"
		
		# atom boxes
		s += "      <g font-size=\"15px\" fill=\"white\" text-rendering=\"optimizeLegibility\" font-family=\"Helvetica\" shape-rendering=\"geometricPrecision\" stroke=\"rgb(0,0,0)\">\n"
		for a in self.graph.atoms:
			s += "        <rect x=\"%d\" y=\"%d\" fill=\"white\" width=\"15\" height=\"11\" stroke=\"none\" />\n" % (a.x+5,a.y-2)
#				if a.symbol == "C":
#					continue
			s += "        <rect x=\"%d\" y=\"%d\" fill=\"white\" width=\"15\" height=\"16\" stroke=\"none\" />\n" % (a.x-7,a.y-8)
		s += "      </g>\n"
		
		# atoms
		s += "      <g font-size=\"15px\" fill=\"rgb(0,0,0)\" text-rendering=\"optimizeLegibility\" font-family=\"Helvetica\" shape-rendering=\"geometricPrecision\" stroke=\"rgb(0,0,0)\">\n"
		for a in self.graph.atoms:
			s += "        <text font-size=\"11px\" fill=\"blue\" xml:space=\"preserve\" x=\"%d\" y=\"%d\" stroke=\"none\">%s</text>\n" % (a.x+6,a.y+8,str(a.id))
#				if a.symbol == "C":
#					continue
				
			if a.symbol == "C":
				color = "black"
			elif a.symbol == "N":
				color = "blue"
			elif a.symbol == "O":
				color = "red"
			elif a.symbol == "P":
				color = "brown"
			elif a.symbol == "S":
				color = "darkgreen"
			else:
				color = "green"
			
			s += "        <text fill=\"%s\" xml:space=\"preserve\" x=\"%d\" y=\"%d\" stroke=\"none\">%s</text>\n" % (color, a.x-5,a.y+5,a.symbol)
		s += "      </g>\n"
		
		# end
		s += "  </g>\n</svg>\n"
		
		f = open(self.ligand + ".svg", "w")
		f.write(s)
		f.close()
	
	
	def __str__(self):
		return self.mol_id
	

