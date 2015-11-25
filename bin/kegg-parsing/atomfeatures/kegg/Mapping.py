# -*- coding: iso-8859-15 -*-

import re, sys
from Reaction import *

class Mapping:
	
	# pit?is menn? niin ett? otetaan reaktio-olio sis??n, sielt? rpairit
	# luetaan sitten rpair-filest? sopivista kohdista...
	
	
	keggmap = "/group/home/icomic/data/kegg/ligand/LATEST/rpair"
	rpairdata = None
	
	def __init__(self, reaction):
		self.reaction = reaction
		self.reaction.mapping = self
		self.mapfile = None
		self.empty = True
		
		self.mappingtype = None
		# Atom -> Atom
		self.mapping = {}
		
		
#		# no rpair given, reading from default location
#		if not rpair:
#			self.mapfile = Mapping.mapfile
#			self.read_kegg_mapping()
#		# rpair filename given, reading from there
#		elif type(rpair) == str:
#			self.mapfile = rpair
#			self.read_kegg_mapping()
#		# rpair data structure given, using it
#		elif type(rpair) == dict:
#			self.read_kegg_mapping_from_struct(rpair)
	
	
#	@staticmethod
	def read_rpair(self, filename = None):
		fn = None
		
		if not filename:
			fn = Mapping.keggmap
		else:
			fn = filename
		
		f = open(fn)
		lines = f.readlines()
		f.close()
		
		# remove everything useless,
		# only leave
		#  - ENTRY
		#  - REACTION
		#  - NAME
		#  - ALIGN
		
#		# reaction -> rpair {R00010 -> [RP00424, RP09422,...] }
#		reaction_to_rpair = {}
		# rpair -> block {RP00400 -> ["ENTRY   RP00400", "REACTION xxxx",...] }
		rpairs = {}
		lastkeyword = None
		entry = None
		for line in lines:
			line = line.rstrip()
			
			if not line.startswith(" "):
				lastkeyword = line.split()[0]
			
			if line.startswith("ENTRY "):
				entry = line.split()[1]
				rpairs[entry] = []
				rpairs[entry].append(line)
			elif line.startswith("NAME"):
				rpairs[entry].append(line)
			elif line.startswith("REACTION") or (lastkeyword == "REACTION" and line.startswith(" ")):
				rpairs[entry].append(line)
			elif line.startswith("ALIGN") or (lastkeyword == "ALIGN" and line.startswith(" ")):
				rpairs[entry].append(line)
		
		Mapping.rpairdata = rpairs
		
		if self.mapping is not None:
			self.empty = False

#		
#		return rpairs
	
	
	def read_kegg_mapping(self, file):
		if not Mapping.rpairdata:
			self.read_rpair(file)
		
		self.mappingtype = "KEGG"
		
		# read the mapping using a given struct
		
		pairs = []
		for rp in self.reaction.rpairs:
#			print rp
			if not Mapping.rpairdata.has_key(rp):
				continue
			
			align = False
			for line in Mapping.rpairdata[rp]:
#				print line
				if line.startswith("NAME"):
					mol1,mol2 = line.split()[1].split("_")
				elif line.startswith("ALIGN"):
					align = True
				elif line.startswith(" ") and align:
					if line.split()[0].strip().isdigit():
						pair = line.split()[1:3]
						a = int(pair[0].split(":")[0])
						b = int(pair[1].split(":")[0])
						pairs.append( (mol1+":"+str(a),mol2+":"+str(b)) )
#						print pairs[-1]
		
#		for x in pairs:
#		print pairs
		
		# RPAIR and COMPOUND/mol files have sometimes different atom descriptors,
		# i.e. the numbering of atoms is different
		# -> check if the numbering doesn't match and abort if doesn't
		#
		
		
		self.construct_mapping(pairs)
	
	
	def complete(self):
		return self.completeness()[0] == self.completeness()[1]
	
	
	def completeness(self):
		mapcount = 0
		total = 0
		for m in self.reaction.subs:
			for a in m.graph.atoms:
				if a.symbol != "H":
					if self.mapping.has_key(a):
						mapcount += 1
					total += 1
		
		return (mapcount, total)
	
	
	
	def construct_mapping(self, pairs):
		# cases:
		# 1) A <-> B
		# 2) A + B <-> C
		# 3) A + A <-> C
		# 4) A <-> B + C
		# 5) A <-> B + B
		#
		# construct mapping by assigning atoms to the pairs
		# problems:
		#  - duplicate compounds share a pairing
		
		mapped = set()
		
		for p in pairs:
			lhsatoms = self.reaction.get_atoms(p[0].split(":")[0],p[0].split(":")[1])
			rhsatoms = self.reaction.get_atoms(p[1].split(":")[0],p[1].split(":")[1])
			
			for lhs in lhsatoms:
				for rhs in rhsatoms:
					if lhs not in mapped and rhs not in mapped:
						self.mapping[lhs] = rhs
						self.mapping[rhs] = lhs
						mapped.add(lhs)
						mapped.add(rhs)
						break
		
		if self.mapping is not None:
			self.empty = False
		
	
	
	
	def read_astar_mapping(self, mapfile = None):
#		print self.reaction.reaction_id
		if not mapfile:
			mapfile = "/fs/home/mqheinon/duuni/atomisota/atom-mappings/astar++/" + self.reaction.reaction_id + ".txt"
		
#		print mapfile
		
		try:
			f = open(mapfile)
		except:
			self.mapping = None
			return
		
		self.mappingtype = "ASTAR"
		
		lines = f.readlines()
		f.close()
		
		# extract first mappings lines from file
		maplines = []
		for line in lines:
			line = line.strip()
			
			if line.startswith("EQUATION"):
				# take the molecules of the reaction
				mols = line.split()[1:]
				mols = [x for x in mols if x.startswith("C")]
			if line.startswith("INDICES"):
				# take indices of the reaction compounds
				indices = line.split()[1:]
				indices = [int(x)-1 for x in indices if x.isdigit()]
			
			if re.search(r'[A-Za-z*?#]+\s+\d+\s+\d+\s+[A-Za-z*?#]+\s+\d+\s+\d+', line):
				maplines.append(line.split())
			if line.startswith("BONDS"):
				break
		
		# reorder the mols by indices
		molmap = [ (indices[i], mols[i]) for i in range(len(mols))]
		molmap.sort()
		mols = [m for i,m in molmap]
		
		# map the mols
		compounds = self.reaction.subs + self.reaction.prods
		for i,m in enumerate(mols):
			for c in compounds:
				if m == c.mol_id and c not in mols:
					mols[i] = c
					break
		
#		print maplines
		
		for line in maplines:
#			print line
			index1 = int(line[2])-1
			index2 = int(line[5])-1
			a = int(line[1])
			b = int(line[4])
			
#			print "ok", a,b,index1,index2
#			print mols[index1].graph
			# we have atom index 'a' of molecule from index 'index1'
			lhs = mols[index1].get_atom(a)
			rhs = mols[index2].get_atom(b)
			self.mapping[lhs] = rhs
		
		if self.mapping is not None:
				self.empty = False
		
		
		


	def read_arita_mapping(self):
		self.mappingtype = "ARITA"
		self.mapping = None
		
		# TODO!!!
		
		print "arita"
		
		f = open(arg)
		lines = f.readlines()
		f.close()
		
		eccode = ""
		formula = ""
		source = ""
		target = ""
		
		mapping = {}
		
		
		for line in lines:
			line = line.strip()
			
			# start
			if "REACTION" in line:
				eccode = ""
				formula = ""
				mapping = {}
			
			# end
			elif "///" in line:
				mappings[eccode] = mapping

#				print eccode
#				for k,v in mappings[eccode].items():
#					print k, v
			
			elif re.match(r'\d+[.][\d|-]+[.][\d|-]+[.][\d|-]+:', line):
				(eccode, formula) = line.split(":")
			
			elif "SOURCE" in line:
				source = line.split()[1]
				
				if not mapping.has_key(source):
					mapping[source] = {}
			
			elif "TARGET" in line:
				target = line.split()[1]
				
				if not mapping.has_key(target):
					mapping[target] = {}
			
			elif ";" in line:
				# map atoms of 'source' to atoms of 'target'
				atompairs = line.split(";")
				atompairs = [x for x in atompairs if x]
				atompairs = [x.strip().split(",") for x in atompairs]
				
#				mapping[molnum][atomnum] = (Cxxxx,3)
				for ap in atompairs:
					mapping[source][ap[0]] = (target, ap[1])
					mapping[target][ap[1]] = (source, ap[0])
				
			
			elif line == "":
				source = ""
				target = ""
		
		if self.mapping is not None:
			self.empty = False
		
	
	
	def compute_cost(self, bc = False):
		# go through all atom pairs (sub_atom, prod_atom)
		#  if bond exists, but not in mapping -> +1
		#  if bond doesn't exist, but in mapping -> +1
		#  if bond exists in both, but of different strength -> +1 (maybe)
		
#		print self
		
		subatoms = []
		for mol in self.reaction.subs:
			for a in mol.graph.atoms:
				if a.symbol != "H":
					subatoms.append(a)
		prodatoms = []
		for mol in self.reaction.prods:
			for a in mol.graph.atoms:
				if a.symbol != "H":
					prodatoms.append(a)
		
#		print "subatoms:", len(subatoms), "prodatoms:", len(prodatoms)
		
		cost = 0
		
		for ind1,a1 in enumerate(subatoms[:-1]):
			for ind2,a2 in enumerate(subatoms[ind1+1:]):
#				print a1, a2
				if bc and a1 in a2.GetAtomNeighbors() and self.mapping[a1] in self.mapping[a2].GetAtomNeighbors():
					origbond = None
					mapbond = None
					
					for b in a1.GetBondNeighbors():
						if b in a2.GetBondNeighbors():
							origbond = b
							break
					
					for b in self.mapping[a1].GetBondNeighbors():
						if b in self.mapping[a2].GetBondNeighbors():
							mapbond = b
							break
					
					if origbond.type != mapbond.type:
						cost += 1
						
				
				elif a1 in a2.GetAtomNeighbors() and self.mapping[a1] not in self.mapping[a2].GetAtomNeighbors():
#					print a1, "is neighbor to:"
#					for ne in a2.GetAtomNeighbors():
#						print ne
#					print a1,a2, "bond broken"
					
					cost += 1
				elif a1 not in a2.GetAtomNeighbors() and self.mapping[a1] in self.mapping[a2].GetAtomNeighbors():
#					print a1,a2, "bond created"
					
					cost += 1
					
		
		return cost
	
	def __str__(self):
		s = ""
		s += self.reaction.equation + "\n"
		for m in self.reaction.subs:
			s += m.mol_id + ":\n"
			for a in m.graph.atoms:
				if self.mapping.has_key(a):
					s += " %s -> %s (%s)\n" % (str(a), str(self.mapping[a]), self.mapping[a].graph.molecule.mol_id)
				else:
					s += " %s ->\n" % (str(a))
		return s.strip()
		
		
		
		
		
		
		


