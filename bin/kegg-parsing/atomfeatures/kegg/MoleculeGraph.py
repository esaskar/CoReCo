# -*- coding: iso-8859-15 -*-

from Atom import *
from Bond import *

class MoleculeGraph:
	SOURCE = 1
	TARGET = 2
	
	ATOMS = ["C","H","N","S","P","O"]
	BONDS = {"1":"-",
        	 "2":"=",
	         "3":"~"}
	
	def __init__(self, molecule = None):
		self.atoms = []
		self.bonds = []
		self.molecule = molecule
		self.compoundatoms = {}
		self.compoundbonds = {}
		self.origatoms = {}
		self.origbonds = {}
	
	
#	# luo uuden irrallisen solmun
#	def AddAtom(self):
#		self.atoms.append( Atom(graph=self) )

	# lisää luodun solmun verkkoon
	def AddAtom(self, atom=None):
		if atom:
			self.atoms.append(atom)
		else:
			self.atoms.append( Atom(symbol="X", graph=self) )

	# luo uuden solmun ja kiinnittää sen kaarella 'to_node':een
	def AddAtomTo(self, to_atom):
		if not to_atom:
			return
		
		newatom = Atom(graph=self) # luodaan uusi solmu
		
		b = Bond(newatom, to_atom) # luodaan kaari 'newnode' -> 'to_node'

		# lisätään molempien solmujen kaarilistaan
		to_atom.AddBond(b)
		newatom.AddBond(b)
		
		# lisätään verkon listoihin uusi solmu ja kaari
		self.atoms.append(newatom)
		self.bonds.append(b)
		
	# Lisätään kaari olemassaolevien solmujen väliin	
	def AddBondTo(self, source, target):
		if not source or not target:
			return
		
		b = Bond(source,target)
		
		# lisätään kaari molempien solmujen kaarilistaan
		source.AddBond(b)
		target.AddBond(b)
		
		# lisätään verkon listaan uusi kaari
		self.bonds.append(b)
	
	# Lisää olemassaolevan kaaren
	def AddBond(self, bond):
		if not bond:
			return
		
		# kiinnitetään kaari päätysolmuihin
		bond.source.AddBond(bond)
		bond.target.AddBond(bond)
	
		self.bonds.append(bond)
	
#	# Lisää olemassa olevan Graph-olion tänne
#	# Kuvaa solmut ja kaaret uusiin id-arvoihin
#	def AddMoleculeGraph(self, other):
#		aidcounter = max( [x.id for x in self.atoms] )
#		bidcounter = max( [x.id for x in self.bonds] )
#		
#		for n in other.atoms:
#			n.id += aidcounter
#		
#		for e in other.bonds:
#			e.id += bidcounter
#		
#		self.atoms += other.atoms
#		self.bonds += other.bonds
	
	
	def __str__(self):
		s = "MoleculeGraph: " + str(len(self.atoms)) + " atoms, " + str(len(self.bonds)) + " bonds\n"
		
		for a in self.atoms:
			s += str(a) + "\n"
		for b in self.bonds:
			s += str(b) + "\n"
		
		return s

	def __len__(self):
		return len(self.atoms)
	
	
	###########################################################################
	###########################################################################
	#
	# mega graph code below
	#
	#
	
	# adds another moleculegraph, with mappings
	def AddMoleculeGraph(self, mg, mapping):
		# we have to take care of:
		# - atoms
		# - bonds
		# - compoundatoms
		# - compoundbonds
		# - origatoms
		# - origbonds
		
		aidcounter = max([x.id for x in self.atoms])+1 if self.atoms else 0 
		bidcounter = max([x.id for x in self.bonds])+1 if self.atoms else 0 
		
		self.compoundatoms[mg.molecule.ligand] = []
		self.compoundbonds[mg.molecule.ligand] = []
		
		# go through atoms on the rhs side
		for a in mg.atoms:
			if a in mapping:
				self.origatoms.setdefault(mapping[a], []).append(a)
				self.compoundatoms.setdefault(mg.molecule.ligand, []).append(a)
			else:
				newatom = Atom(symbol=a.symbol, id=aidcounter, graph=self)
				self.AddAtom(newatom)
				aidcounter += 1
				mapping[a] = newatom
				self.origatoms.setdefault(newatom, []).append(a)
				self.compoundatoms.setdefault(mg.molecule.ligand,[]).append(newatom)
		
		for b in mg.bonds:
			newsrc = mapping[b.source]
			newtgt = mapping[b.target]
			mbond = newsrc.GetBond(newtgt)
			
			if mbond:
				self.origbonds.setdefault(mbond,[]).append(b)
				self.compoundbonds.setdefault(mg.molecule.ligand,[]).append(mbond)
				mapping[b] = mbond
			else:
				newbond = Bond(newsrc, newtgt, b.type, bidcounter, self)
				self.AddBond(newbond)
				bidcounter += 1
				mapping[b] = newbond
				self.origbonds.setdefault(newbond,[]).append(b)
				self.compoundbonds.setdefault(mg.molecule.ligand,[]).append(newbond)
	
	
	def output(self, name=None):
		if not name:
			name = self.molecule.ligand + ".mg"
		
		s = ""
		
		for a in self.atoms:
			coords = {}
			for origatom in self.origatoms[a]:
				c = origatom.graph.molecule.ligand
				coords[c] = (origatom.x,origatom.y)
			
			coordset = ["%s:%.3f|%.3f" % (c,x,y) for c,(x,y) in coords.items()]
			
			s += "ATOM %d %s %s\n" % (a.id, a.symbol, ",".join(coordset) )
		
		for b in self.bonds:
			s += "BOND %d %d %d %d\n" % (b.id, b.source.id, b.target.id, int(b.type))
		
		for c in self.compoundatoms:
			atoms = self.compoundatoms[c]
			bonds = self.compoundbonds[c]
			s += "COMPOUND %s [%s] [%s]\n" % (c, ",".join([str(a.id) for a in atoms]), ",".join([str(b.id) for b in bonds]))
		
		f = open(name, "w")
		f.write(s)
		f.close()
		
	
	def modular_product_graph(self, other):
		Atom.idcounter = 0  # reset counters to get zero-based id's
		Bond.idcounter = 0
		
		mpg = MoleculeGraph()
		pairs = [(a1,a2) for a1 in self.atoms for a2 in other.atoms if a1.symbol == a2.symbol]
		for a1,a2 in pairs:
			mpg.AddAtom()
			mpg.atoms[-1].a1 = a1
			mpg.atoms[-1].a2 = a2
		
		
#		print "modular graph"
#		for a in mpg.atoms:
#			print a, a.a1, a.a2
		print "found", len(mpg.atoms)
		
		for v in mpg.atoms:
			for u in mpg.atoms:
				if v.a1.IsNeighbor(u.a1) == v.a2.IsNeighbor(u.a2):
					mpg.AddBondTo(u,v)
#					print "added bond", u.a1, u.a2, "<->", v.a1, v.a2
#				elif not v.a1.IsNeighbor(u.a1) and not v.a2.IsNeighbor(u.a2):
#					mpg.AddBondTo(u,v)
#					print "added bond", u.a1, u.a2, "<->", v.a1, v.a2
		
		print "graph1", self
		for a in self.atoms:
			print a
		for b in self.bonds:
			print b
		
		print
		print "graph2", other
		for a in other.atoms:
			print a
		for b in other.bonds:
			print b
		
		
		print 
		print "modular graph"
		for a in mpg.atoms:
			print a, a1, a2
#		for b in mpg.bonds:
#			print b
		
		
		return mpg
	





















