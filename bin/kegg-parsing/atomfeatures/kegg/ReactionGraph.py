# -*- coding: iso-8859-15 -*-
#
#
# Reaction Graph is a supergraph containing information about reactions and compounds
# There each reactiongraph has a parent "MetabolicNetwork" which contains a bunch of ReactionGraphs,
# as well as a set of original reactions and compounds
# 
#
# Atoms have no extra fields 
# Bonds have:
# - reactions { reaction-str -> [+1|0|-1] }, ie. "R00882":+1 means a new bond
# 
# This class has:
# - atoms [obj]
# - bonds [obj]
# - reaction-obj (original reaction)
# - compounds { compound-str -> [newatoms] }, ie. "C00083": [obj,ojb,obj...]
# - reactions { reaction-str -> [newatoms] }
# - origatoms { newatom -> oldatom } (pointer to original atoms)
# - origbonds { newbond -> oldbond }
#
#

import math, time
import kegg
from Atom import *
from Bond import *
from kegg import VF

class Point:
	pass


class ReactionGraph:
	SOURCE = 1
	TARGET = 2
	ATOMS = ["C","H","N","S","P","O"]
	BONDS = {"1":"-", "2":"=", "3":"~"}
	
	def __init__(self, reaction = None):
		self.id = 0
		self.atoms = [] # set of new atoms for rgraph
		self.bonds = [] # set of new bonds
		self.reaction = reaction
		self.compoundatoms = {} # points to the new atoms inside the mrg, { mol -> [[newatoms],[newatoms],..] }
		self.compoundbonds = {} # points to the new bonds inside the mrg, { mol -> [[newbonds],[newbonds],..] }
		self.reactions = {} # points to the new atoms inside the mrg, { reac -> [newatoms] }
		self.origatoms = {} # points to the actual atoms from inside the mrg, { newatom -> [origatoms] }
		self.origbonds = {} # points to the actual bonds from inside the mrg, { newbond -> [origbonds] }
	
	
	# produces a description of the reaction graph
	def info(self):
		s = "Reaction graph %s (%d atoms, %d bonds)\n" % (self.reaction.ligand if self.reaction else "", len(self.atoms), len(self.bonds))
		
		for a in self.atoms:
			mols = []
			for mol,pieces in self.compoundatoms.items():
				for atoms in pieces:
					if a in atoms:
						mols.append(mol)
			reacs = []
			for reac,atoms in self.reactions.items():
				if a in atoms:
					reacs.append(reac)
			
			s += "%s -> %s { comps: %s } { reacs: %s }\n" % (str(a), ",".join([str(ne) for ne in a.GetAtomNeighbors()]), ",".join(mols), ",".join(reacs))
		
		for b in self.bonds:
			mols = []
			for mol,pieces in self.compoundbonds.items():
				for bonds in pieces:
					if b in bonds:
						mols.append(mol)
			
			s += "%s -> { %s } { comps: %s }\n" % (str(b), ",".join([(str(r.ligand) if isinstance(r,kegg.Reaction) else r) + ":" + str(v) for r,v in b.reactions.items() if v == 1 or v == -1]), ",".join(mols)) 
		
		for mol,pieces in self.compoundatoms.items():
			s += "%s: " % (mol)
			for atoms in pieces:
				s += "["
				s += ",".join([str(a.id) for a in sorted(atoms)])
				s += "]"
			s += "\n"
#			s += "%s: %s\n" % (mol, ",".join([str(a.id) for a in sorted(ats)]))
		for r,ats in self.reactions.items():
			s += "%s: [ %s ]\n" % (r, ",".join([str(a.id) for a in sorted(ats)]))
		
		return s
	
	
	def read(self, mapfile = None):
#		print self.reaction.ligand
		# read mappings
		
		self.reaction.read_mapping("ASTAR", mapfile)
		if not self.reaction.mappings[0].mapping:
#			print "no mappings"
			self = None
			return
		
		# get sub/prod atoms, mapping both ways, compound list
		subatoms = [a for sub in self.reaction.subs if sub.graph for a in sub.graph.atoms ]
		prodatoms = [a for prod in self.reaction.prods if prod.graph for a in prod.graph.atoms]
		mapping = self.reaction.mappings[-1].mapping
		revmapping = dict( [(rhs,lhs) for lhs,rhs in mapping.items()] )
		newatoms = {}
		
		aids = 1 # atom id counter
		bids = 1
		
#		for k,v in mapping.items():
#			print k,v
#		print
#		for k,v in revmapping.items():
#			print k,v
		
		
		# go through subs and prods and create new atoms based on these
		for a in subatoms:
			newatom = Atom(a.symbol, id=aids, graph=self)
			aids += 1
			newatom.x = a.x
			newatom.y = a.y
			self.AddAtom(newatom)
			newatoms[a] = newatom
		
		#
		# the self.compoundatoms and self.compoundbonds contain { "C00031":[[atoms],[atoms],...] }
		# if there are multiple same reactant (eg. 2 C00031 <=> C00634 + C03259), we need to map both
		# C00031's into indices in compoundatoms-list of atoms
		# this id done with 'molindex{}'
		#
		# TODO: change .compoundatoms and .compoundbonds to have [[atoms],[atoms]...] structure!!!
		# here change so that we have multiple compounds when eg. 2 C00031 <=> C00634 + C03259
		#
		
		molindex = {}
		
		for a in subatoms:
			mol = a.graph.molecule
			reac = a.graph.molecule.reaction
			newatom = newatoms[a]
			# here for "2 C00031 <=> ******", put "C00031":[[atoms],[atoms]]
			self.compoundatoms.setdefault(mol.ligand, [])
			self.compoundbonds.setdefault(mol.ligand, [])
			self.reactions.setdefault(reac.ligand, [])
			
			# if mol is not seen before, add new atomlist and save its index in the outer list
			if mol not in molindex:
				molindex[mol] = len(self.compoundatoms[mol.ligand])
				self.compoundatoms[mol.ligand].append( [newatom] )
				self.compoundbonds[mol.ligand].append( [] )
			# mol is seen before, simply append
			else:
				index = molindex[mol]
				self.compoundatoms[mol.ligand][index].append(newatom)
			
			self.origatoms.setdefault(newatom, []).append( a )
			
			if newatom not in self.reactions[reac.ligand]:
				self.reactions[reac.ligand].append( newatom )
		
		for a in prodatoms:
#			print a
			mol = a.graph.molecule
			newatom = newatoms[revmapping[a]]
			self.compoundatoms.setdefault(mol.ligand, [])
			self.compoundbonds.setdefault(mol.ligand, [])
			
			# if mol is not seen before, add new atomlist and save its index in the outer list
			if mol not in molindex:
				molindex[mol] = len(self.compoundatoms[mol.ligand])
				self.compoundatoms[mol.ligand].append( [newatom] )
				self.compoundbonds[mol.ligand].append( [] )
			# mol is seen before, simply append
			else:
				index = molindex[mol]
				self.compoundatoms[mol.ligand][index].append(newatom)
			
			self.origatoms.setdefault(newatom, []).append( a )
		
		
		# go through pairs of atoms, add bonds to reaction graph
		for ind1,a1 in enumerate(subatoms[:-1]):
			for a2 in subatoms[ind1+1:]:
				atom1mol = a1.graph.molecule.ligand
				mappedatom1mol = mapping[a1].graph.molecule.ligand
				
				index = molindex[a1.graph.molecule]
				mappedindex = molindex[mapping[a1].graph.molecule]
				
				# intact bond
				if a1 in a2.GetAtomNeighbors() and mapping[a1] in mapping[a2].GetAtomNeighbors():
					b = a1.GetBond(a2)
					
					source = newatoms[b.source]
					target = newatoms[b.target]
					
					newbond = Bond(source, target, b.type, bids, self)
					bids += 1
					newbond.reactions = { self.reaction : 0 } # set a new field
					
					self.AddBond( newbond )
					self.origbonds.setdefault(newbond, []).append(b)
					self.compoundbonds[atom1mol][index].append( newbond )
					self.compoundbonds[mappedatom1mol][mappedindex].append( newbond )
				
				# cleaved bond
				elif a1 in a2.GetAtomNeighbors() and mapping[a1] not in mapping[a2].GetAtomNeighbors():
					b = a1.GetBond(a2)
					
					source = newatoms[b.source]
					target = newatoms[b.target]
					
					newbond = Bond(source, target, b.type, bids, self)
					bids += 1
					newbond.reactions = { self.reaction : -1 } # set a new field
					
					self.AddBond(newbond)
					self.origbonds.setdefault(newbond, []).append(b)
					self.compoundbonds[atom1mol][index].append( newbond )
					
				
				# new bond
				elif a1 not in a2.GetAtomNeighbors() and mapping[a1] in mapping[a2].GetAtomNeighbors():
					b = mapping[a1].GetBond(mapping[a2])
					
					source = newatoms[a1]
					target = newatoms[a2]
					
					newbond = Bond(source, target, b.type, bids, self)
					bids += 1
					newbond.reactions = { self.reaction : 1 } # set a new field 
					
					self.AddBond(newbond)
					self.origbonds.setdefault(newbond, []).append(b)
					self.compoundbonds[mappedatom1mol][mappedindex].append( newbond )
	
	# add a node atom to the graph
	def AddAtom(self, atom):
		self.atoms.append(atom)
	
	# adds a new atom and attaches it to 'to_atom'
	def AddAtomTo(self, to_atom):
		if not to_atom:
			return
		
		newatom = Atom(graph=self) # new atom
		b = Bond(newatom, to_atom) # new bond
		
		# add to both bonds edge lists
		to_atom.AddBond(b)
		newatom.AddBond(b)
		
		self.atoms.append(newatom)
		self.bonds.append(b)
	
	# adds existing bond
	def AddBond(self, bond):
		if not bond:
			return
		
		# attaches the bond to end points
		bond.source.AddBond(bond)
		bond.target.AddBond(bond)
	
		self.bonds.append(bond)
	
	
	
	
	# (1) compute all automorphisms between shared reactant atoms
	#
	# use simple isomorphism 
	#
	def automorphisms(self, atoms1, atoms2):
		# atoms1 = selfatoms of reactant
		# atoms2 = rgatoms of reactant
		pass
	
	
	# (2) extend the automorphism mapping with all existing matching pairs, the automorphism is fixed
	#
	# state space representation algorithm:
	# state = (map, bordercands), i.e. ( {1:1, 2:2, 3:3}, [(4,4), (4,5), (5,4), (5,5), (8,7)] )
	#
	# use dfs recursive function to:
	# - choose one candidate pair (a,b)
	# - add (a,b) to map
	# - remove from cands all (a,*) and (*,b) pairs
	# - add to cands all neighborpairs of (a,b) which are NOT in map already, 
	#    AND which don't contradict automorphism 
	#    (i.e. if (2,3) in automorphism, only accept (2,3), (4,4), etc...
	#     not (2,4), (5,3), etc... )
	#
	# Ie.: additionally if cands contains automorphism pair (a,b), remove again (a,*) and (*,b)
	# (always choose the direction of automorphism)
	#
	# recursion returns when border is empty
	# a data structure saves only the results with maximum length (i.e. leafs of dfs tree of maximal depth)
	#
	# example: 
	# here next round could be: (a,b) = (4,4)
	# new cands is [(5,5), (8,7)]        # after removal of overlapping cands
	#              [(5,5), (8,7), (6,6), (6,9), (9,6), (9,9)]  # after addition of new neighpairs
	#              [(5,5), (8,7), (6,6), (9,9)]  # after removal of automorphism pairs, if (6,6) was in automorph
	# new map is {1:1, 2:2, 3:3, 4:4}
	#
	#
	def extendmappings(self, automaps):
		self.resultsize = 0
		self.results = []
		self.iters = 0
		
		for automap in automaps:
			rhsborder = [ne for a in automap.keys() for ne in a.GetAtomNeighbors() if ne not in automap]
			lhsborder = [ne for a in automap.values() for ne in a.GetAtomNeighbors() if ne not in automap.values()]
			
			self.dfs(automap, lhsborder, rhsborder)
		
		print self.iters,
		
		return self.results # extended mappings
	
	
	def dfs(self, map, lhs, rhs):
		self.iters += 1
#		prefix = " " * (len(map)-28)
		
		# (1) first loop over atom pairs and try to deepen recursion
		# (2) then check if we couldn't recurse -> result?
		
		recurse = False
		bfirstlist = [rhs[0]] if rhs else []
		for a,b in [(x,y) for x in lhs for y in bfirstlist]: # take all a's and only first b, don't go at all if b empty
			# atom criteria
			if a.symbol != b.symbol:
				continue
			
			# connectedness-criteria: for pair (a,b), we have M[neighs(a)] = neighs(b)
			# i.e. the neighbors of a have to be mapped to neighbors of b
			# 
			# change the criteria: only one neighbor of a has to be mapped to a neighbor of b, not all
			a_mapped_neighs = [ ne for ne in a.GetAtomNeighbors() if ne in map.values() ]
			b_mapped_neighs = [ ne for ne in b.GetAtomNeighbors() if ne in map.keys() ]
			
			for x in b_mapped_neighs:
				if map[x] in a_mapped_neighs: # success!
					break
			else:
				continue
			
			recurse = True
			
			newmap = map.copy()
			newmap[b] = a          # add to map
			newlhs = [x for x in lhs if x != a]
			newrhs = [x for x in rhs if x != b]
			
			a_neighs = [ne for ne in a.GetAtomNeighbors() if ne not in map.values() and ne not in lhs ]  # neighs not mapped yet
			b_neighs = [ne for ne in b.GetAtomNeighbors() if ne not in map.keys()   and ne not in rhs ]
			
			newlhs += a_neighs
			newrhs += b_neighs
			
			self.dfs(newmap, newlhs, newrhs)
		
		if not recurse:
			if len(map) > self.resultsize:   # clear results list
				self.results = [map]
				self.resultsize = len(map)
#			elif len(map) == self.resultsize:   # add new result
#				self.results.append(map)
	
	
	# (3) complete the extended automorphism mapping with regions of rg (nonmatching parts)
	#
	# fill with nones
	#
	def completemappings(self, rg, extendedmaps):
		for m in extendedmaps:
			for a in rg.atoms:
				if a not in m:
					m[a] = None
		
		return extendedmaps
	
	
	# generates optimal mapping from rg to self (rg -> self)
	def optimal_atom_mappings(self, rg):
		# reactions: A + B <=> C + D
		#            E + B <=> G + H
		#            A + E <=> C + F
		# 
		# first and second have a single reactant in common: map the B's of both reactions together
		#  and expand the mapping from there
		# first and third have 1 and 1 reactants in common: map A's together, C's together
		#  and expand the mapping to cover also E and F
		# there shouldn't be overlaps 
		#
		
		t = time.time()
		
		shared_reactants = [mol for mol in self.compoundatoms if mol in rg.compoundatoms]
		shared_subs = [mol for mol in shared_reactants if mol in [m.ligand for m in rg.reaction.subs]]
		shared_prods = [mol for mol in shared_reactants if mol in [m.ligand for m in rg.reaction.prods]]
		if len(shared_subs) > len(shared_prods):
			shared_reactants = shared_subs
		elif len(shared_subs) == len(shared_prods):
			if len([a for m in shared_subs for a in self.compoundatoms[m]]) >= len([a for m in shared_prods for a in self.compoundatoms[m]]):
				shared_reactants = shared_subs
			else:
				shared_reactants = shared_prods
		else:
			shared_reactants = shared_prods
		
		print "shared reactant(s):", ",".join(shared_reactants),
		
		print time.time() - t,
		t = time.time()
		
		
		# with 2 shared reactants, where both are twice in the compounds:
		# shared_reactants = ["C00031", "C00001"]
		# sr = "C00031"
		# self.compoundatoms[sr] = [[atoms],[atoms]]
		# rg.compoundatoms[sr] =   [[atoms],[atoms]]
		#
		# the mappings of several shared reactants should complement each other
		# 
		
		compoundmaps = {}
		for sr in shared_reactants:
#			print sr
			compoundmaps[sr] = []
			for i in range(len(self.compoundatoms[sr])):
				for j in range(len(rg.compoundatoms[sr])):
#					print i,j
#					print "selfatoms", map(str, self.compoundatoms[sr][i])
#					print "selfbonds", map(str, self.compoundbonds[sr][i])
#					print "rgatoms", map(str, rg.compoundatoms[sr][j])
#					print "rgbonds", map(str, rg.compoundbonds[sr][j])
					selfatoms = self.compoundatoms[sr][i]
					selfbonds = self.compoundbonds[sr][i]
					rgatoms   = rg.compoundatoms[sr][j]
					rgbonds   = rg.compoundbonds[sr][j]
					
					vf = VF.SubRegionVF( rgatoms, rgbonds, selfatoms, selfbonds )
#					if sr == "C00446":
#						vf.verbose = True
					automaps = vf.isomorph()
#					if sr == "C00446":
#						print len(automaps), automaps
					for m in automaps:
						compoundmaps[sr].append(m)
#					if sr == "C00446":
#						print compoundmaps
		
		# R01900: compoundmaps contains 16 C417 mappings, 1 C1->22 mapping and 1 C1->161 mapping
		#         all in all [ [{},{},{},...], [{}], [{}] ]
		# what we want: { C417:[{},{},...],   C1:[{},{}] }
		
		mappings = []
		stack = [({},0)] # empty map initially
		while stack:
			(m1,d) = stack.pop()
			
#			if rg.reaction.ligand == "R00955":
#				print map(str, m1.values())
			
			mappings.append(m1)
			
			if d == len(shared_reactants): # leaf node
				continue
			
			sr = shared_reactants[d]
			
#			if rg.reaction.ligand == "R00955":
#				print sr
			
			for m2 in compoundmaps[sr]:
				# check for inconsistencies
				# TODO: check that keys AND values don't overlap with each other, with sets perhaps??
				
#				if rg.reaction.ligand == "R00955":
#					print map(str, m2.values())
				
				if len(set(m2.values()) & set(m1.values())) == 0: # intersection empty
					nextm = m1.copy()
					nextm.update(m2)
					stack.append( (nextm,d+1) )
#				else:
#					if rg.reaction.ligand == "R00955":
#						print "ZZZZZZZZZZZ"
#						print map(str, m2.values())
#						print map(str, m1.values())
#						print map(str, set(m2.values()) & set(m1.values()))
		
		# take only longest mappings
		mappings2 = []
		maxlen = 0
		for m in mappings:
			if len(m) > maxlen:
				mappings2 = [m]
				maxlen = len(m)
			elif len(m) == maxlen:
				mappings2.append(m)
		mappings = mappings2
		
		# EXPERIMENTAL!!
		mappings = [mappings[0]]
		
		
		print "found", len(mappings), "automorphisms",
		
		print time.time() - t,
		t = time.time()
		
		mappings = self.extendmappings(mappings)
		mappings = self.completemappings(rg, mappings)
		
		print len(mappings), "mappings",
		
		print time.time() - t,
		
		return mappings
	
	
	# a method to _try_ to match rg to self to see the resulting size,
	# doesn't actually change self
	def MergeSize(self, mappings):
		return max( [len(m) for m in [[(a,b) for a,b in m.items() if a and b] for m in mappings ]] )
	
	
	# we need to get the common reactants of two reaction graphs
	# there should be a single unique compound in: 
	# commoncomp = [mol for mol in self.compounds if mol in rg.compounds][0]
	# then we get the atoms of both sides:
	# selfatoms = self.compounds[commoncomp]
	# rgatoms = rg.compounds[commoncomp]
	# we need a map between, automorphisms cause symmetries to show
	# for map in mappings:
	#   map[selfatoms] -> rgatoms
	# now for each possible automorphism, merge the shared reactant atoms and bonds,
	# and continue using recursive bfs() algorithm mapping the atoms from rhs to lhs
	#
	#
	# Fuses the 'rg' reaction graph into self
	def AddReactionGraph(self, rg, mappings = None):
		# self = current reaction graph
		# rg = reaction graph to be merged
		#
		# for the mapped regions, add to those atoms the mapped atom's molecules and reactions
		# the remaining unmapped regions atoms should be transferred to this graph
		# the mapped bonds should be added to mcs parts (reactions)
		# the unmapped bonds should be transferred
		# for svg the coordinates have to be calibrated,
		#  basically just overlay the two reaction graphs on top of each other
		#  causes collisions, need aritas layout algorithm?
		#
		
		if not mappings:
			mappings = self.optimal_atom_mappings(rg)
		
		aidcounter = max( [x.id for x in self.atoms] ) + 1
		bidcounter = max( [x.id for x in self.bonds] ) + 1
		
		newidstart = aidcounter
		
		# optimal mappings, basically iso/automorphic
		mapping = mappings[0]
		
		# we need to handle:
		# - reactions
		# - compounds
		# - origatoms
		# - origbonds
		# - bond.reactions
		
		# first add atoms not yet on 'self'
		for a in rg.atoms:
			if mapping[a] is None:
				newatom = Atom(symbol=a.symbol, id=aidcounter, graph=self)
				mapping[a] = newatom
				self.AddAtom(newatom)
				aidcounter += 1
		
		# next add bonds 
		# now we have all atoms mapped from rg to self
		# we have created all new atoms
		# we have all atom information updated
		# next we go through all bonds in rg and map them to self and create all new bonds
		#
		for b in rg.bonds:
			src = b.source
			tgt = b.target
			msrc = mapping[src] # can be none
			mtgt = mapping[tgt]
			mb = msrc.GetBond(mtgt)
			
			# bond exists
			if msrc and mtgt and mb:
				mb.reactions.update(b.reactions)
				mapping[b] = mb
			# new bond
			else:
				newbond = Bond(msrc, mtgt, b.type, bidcounter, self)
				newbond.reactions = b.reactions
				# if the new bond attaches new atom to old atom, set the bond as +1
#				if (newbond.source.id < newidstart) != (newbond.target.id < newidstart):
#					newbond.reactions[rg.reaction] = +1
				if newbond.source.id < newidstart or newbond.target.id < newidstart:
					newbond.reactions[rg.reaction] = +1
				
				bidcounter += 1
				self.AddBond(newbond)
				mapping[b] = newbond
		
		# next update information on the atoms
		for rhs in rg.atoms:
			lhs = mapping[rhs]
			
			# origatoms ( newatom -> [origatoms] )
			try:
				self.origatoms[lhs] += rg.origatoms[rhs]
			except:
				self.origatoms[lhs] = rg.origatoms[rhs]
		
		# finally update information on all rg bonds mapped to self
		# origbonds ( newbond -> [origbonds] )
		for rhs in rg.bonds:
			lhs = mapping[rhs]
			
			try:
				self.origbonds[lhs] += rg.origbonds[rhs]
			except:
				self.origbonds[lhs] = rg.origbonds[rhs]
		
		# TODO: a compound can exist in numerous places, the compoundatoms should be { mol: [[atoms],[atoms],..] }
		for c,pieces in rg.compoundatoms.items():
			for atoms in pieces:
				if c in self.compoundatoms:
					indices = {}
					for i in range(len(self.compoundatoms[c])):
						Alist = self.compoundatoms[c][i]
						for a in atoms:
							if mapping[a] in Alist:
								indices[a] = i
					components = len(set(indices.values()))
				
				if c not in self.compoundatoms or components != 1:
					self.compoundatoms.setdefault(c,[]).append([])
					for a in atoms:
						self.compoundatoms[c][-1].append( mapping[a] )
		
		for c,pieces in rg.compoundbonds.items():
			for bonds in pieces:
				if c in self.compoundbonds:
#					matching_areas = list(dict([ (B,b) for B in self.compoundbonds[c] for b in bonds if mapping[b] in B]))
					indices = {}
					for i in range(len(self.compoundbonds[c])):
						Blist = self.compoundbonds[c][i]
						for b in bonds:
							if mapping[b] in Blist:
								indices[b] = i
					components = len(set(indices.values()))
				
				if c not in self.compoundbonds or components != 1:
					self.compoundbonds.setdefault(c,[]).append([])
					for b in bonds:
						self.compoundbonds[c][-1].append( mapping[b] )
		
		for r,ats in rg.reactions.items():
			if r not in self.reactions:
				self.reactions[r] = []
				for a in ats:
					self.reactions[r].append( mapping[a] )
	
	
	
	
	def graphviz(self, fname):
		s = "graph g {\n\n"
		
		for a in self.atoms:
			if a.symbol == "C":
				color = "black"
			elif a.symbol == "O":
				color = "red"
			elif a.symbol == "N":
				color = "blue"
			elif a.symbol == "S":
				color = "brown"
			elif a.symbol == "P":
				color = "yellow"
			else:
				color = "green"
			s += '"%s %d" [color="%s"];\n' % (a.symbol, a.id, color)
		s += "\n"
		
		for b in self.bonds:
			comps = len([c for c,bondsets in self.compoundbonds.items() if b in bondsets[0]])
			length = 1.0
			if comps > 0:
				length = math.log(comps) + 1.0
			s += '"%s %d" -- "%s %d" [len=%.2f];\n' % (b.source.symbol, b.source.id, b.target.symbol, b.target.id, length)
		s += "\n"
		s += "}\n"
		
		f = open(fname, "w")
		f.write(s)
		f.close()
	
	def output_molfile(self, fname):
		# basic format:
		# three header lines (can be empty)
		# stats-line, with format AAABBB  0  0  0  0            999 V2000
		# atom-line,  with format  +XXX.XXXX +YYY.YYYY +ZZZ.ZZZZ S   0  0  0  0  0  0  0  0  0  0  0  0
		# bond-line,  with format IIIJJJBBB  0  0  0  0
		# endline,    with format M  END
		
		atomcount = len(self.atoms)
		bondcount = len(self.bonds)
		
		try:
			ec = " ".join(self.reaction.ec)
		except:
			ec = ""
		
		s = "%s %s\n\n\n" % (self.reaction.ligand, ec) # header
		s += "%s%s  0  0  0  0            999 V2000\n" % (str(atomcount).rjust(3), str(bondcount).rjust(3))
		for a in self.atoms:
			s += "    0.0000    0.0000    0.0000 %s   0  0  0  0  0  0  0  0  0  0  0  0\n" % (a.symbol)
		for b in self.bonds:
			rtype =  str(b.reactions[self.reaction])
			s += "%s%s%s%s  0  0  0\n" % (str(b.source.id).rjust(3), str(b.target.id).rjust(3), str(b.type).rjust(3), rtype.rjust(3))
		s += "M  END\n"

		f = open(fname, "w")
		f.write(s)
		f.close()
	
	def output(self, fname):
		# outputs the reaction graph into a file
		# 
		# the format
		#
		# ATOM ID SYMBOL X.XX Y.YY
		# ...
		# BOND ID SOURCE TARGET TYPE
		# ...
		# COMPOUND ID [AIDS] [BIDS]
		# ...
		# REACTION ID [AIDS] [BIDS]
		# ...
		# REACTIONCHANGE BID RID:{+1|-1}
		# ...
		#
		# example:
		# ATOM 1 C   C00031:X.XX|Y.YY
		# ATOM 2 O 
		# ATOM 3 N
		# BOND 1 1 2 1
		# COMPOUND C00031 [2,3,4] [2,3]
		# REACTION R00421 [3,4] [5,2]
		# REACTIONCHANGE 2 R00021:-1 
		# REACTIONCHANGE 2 R00047:+1 
		# 
		
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
		for c,atomsets in self.compoundatoms.items():
			for i in range(len(atomsets)):
				atoms = self.compoundatoms[c][i]
				bonds = self.compoundbonds[c][i]
				s += "COMPOUND %s [%s] [%s]\n" % (c, ",".join([str(a.id) for a in atoms]), ",".join([str(b.id) for b in bonds]))
		for r,atoms in self.reactions.items():
			bond_ids = []
			for b in self.bonds:
				if r in [reac.ligand for reac in b.reactions]:
					bond_ids.append(str(b.id))
			s += "REACTION %s [%s] [%s]\n" % (r, ",".join([str(a.id) for a in atoms]), ",".join(bond_ids))
		for b in self.bonds:
			for k,v in b.reactions.items():
				if v != 0:
					s += "REACTIONCHANGE %d %s:%s\n" % (b.id, k, str(v))
		
		f = open(fname, "w")
		f.write(s)
		f.close()
	
	def input(self, fname):
		f = open(fname)
		lines = map(str.strip, f.readlines())
		f.close()
		
		# set
		# - self.atoms
		# - self.bonds
		# - self.reactions
		# - self.compoundatoms
		# - self.compoundbonds
		
		for line in lines:
			data = line.split()[1:]
			
#			print line
			
			if line.startswith("ATOM"):
				i = int(data[0])
				s = data[1]
				coords = {}
				
				for xystr in data[2].split(","):
					c = xystr.split(":")[0]
					x,y = xystr.split(":")[1].split("|")
					coords[c] = Point()
					coords[c].x = float(x)
					coords[c].y = float(y)
				
				a = Atom(id=i,symbol=s)
				a.coords = coords
				
				self.AddAtom(a)
				
			elif line.startswith("BOND"):
				i = int(data[0])
				srcid = int(data[1])
				tgtid = int(data[2])
				t = int(data[3])
				source = self.atoms[srcid-1]
				target = self.atoms[tgtid-1]
				b = Bond(id=i,source=source,target=target,type=t)
				b.reactions = {}
				self.AddBond(b)
				
			elif line.startswith("COMPOUND"):
				c = data[0]
				atoms = data[1][1:-1]
				bonds = data[2][1:-1]
				
				# first compoundptr
				if c not in self.compoundatoms:
					self.compoundatoms[c] = [[]]
					self.compoundbonds[c] = [[]]
					for aid in atoms.split(","):
						self.compoundatoms[c][-1].append( self.atoms[int(aid)-1])
					for bid in bonds.split(","):
						if bid:
							self.compoundbonds[c][-1].append(self.bonds[int(bid)-1])
				else:
					self.compoundatoms[c].append([])
					self.compoundbonds[c].append([])
					for aid in atoms.split(","):
						self.compoundatoms[c][-1].append( self.atoms[int(aid)-1])
					for bid in bonds.split(","):
						if bid:
							self.compoundbonds[c][-1].append(self.bonds[int(bid)-1])
			
			elif line.startswith("REACTIONCHANGE"):
				bid = int(data[0])
				ligand = data[1].split(":")[0]
				type = int(data[1].split(":")[1])
				
				self.bonds[bid-1].reactions[ligand] = type
			
			elif line.startswith("REACTION"):
				ligand = data[0]
				atoms = data[1][1:-1]
				bonds = data[2][1:-1]
				
				self.reactions[ligand] = []
				for aid in atoms.split(","):
					self.reactions[ligand].append(self.atoms[int(aid)-1])
				
				for bid in bonds.split(","):
					bid = int(bid)
					self.bonds[bid-1].reactions[ligand] = 0
			
	
	
	
	def compute_borders(self):
		width = max([a.x for a in self.atoms]) - min([a.x for a in self.atoms])
		height = max([a.y for a in self.atoms]) - min([a.y for a in self.atoms])
		
		return width + 200, height + 200
	
	
	# scales the coordinates such that minimum distance between atoms is '24' (default)
	# positions the coords such that left-top corner is (20,20) and right-bottom corner (n,m)
	# also adjusts the all coordinates by 'dx' and 'dy'
	def scale_atoms(self, atomdist=24, dx=0, dy=0):
		def atomdistance(a,b):
			return math.sqrt((a.x-b.x)**2 + (a.y-b.y)**2)
		
		# scale coordinates such that distance between closest atoms is exactly 24px (default)
		left = top = 20
		# find min-distance between any two atoms
		mindist = min([atomdistance(a1,a2) for a1 in self.atoms for a2 in self.atoms if a1 != a2])
		# scale
		for a in self.atoms:
			a.x *= atomdist / mindist
			a.y *= atomdist / mindist
		
		# make all coordinates positive and left-top corner -most atom into (20,20) coordinate
		
		min_x = min([a.x for a in self.atoms])
		max_x = max([a.x for a in self.atoms])
		min_y = min([a.y for a in self.atoms])
		max_y = max([a.y for a in self.atoms])
		
		for a in self.atoms:
			a.x -= (min_x - left - dx)
			a.y -= (min_y - top - dy)
	
	def center_atoms(self):
		left = top = 50
		# make all coordinates positive and left-top corner -most atom into (50,50) coordinate
		min_x = min([a.x for a in self.atoms])
		max_x = max([a.x for a in self.atoms])
		min_y = min([a.y for a in self.atoms])
		max_y = max([a.y for a in self.atoms])
		
		for a in self.atoms:
			a.x -= (min_x - left)
			a.y -= (min_y - top)
		
	
	def scale_compound(self, c, atomdist=24, dx=0, dy=0):
		def atomdistance(a1,a2):
			return math.sqrt((a1.coords[c].x-a2.coords[c].x)**2 + (a1.coords[c].y-a2.coords[c].y)**2)
		
		left = top = 20
		atoms = self.compoundatoms[c][0]
		mindist = min([atomdistance(a1,a2) for a1 in atoms for a2 in atoms if a1 is not a2]) if len(atoms) > 1 else atomdist
		
		for a in atoms:
			a.coords[c].x *= atomdist / mindist
			a.coords[c].y *= atomdist / mindist
		
		min_x = min([a.coords[c].x for a in atoms])
		max_x = max([a.coords[c].x for a in atoms])
		min_y = min([a.coords[c].y for a in atoms])
		max_y = max([a.coords[c].y for a in atoms])
		
		for a in atoms:
			a.coords[c].x -= (min_x - left - dx)
			a.coords[c].y -= (min_y - top - dy)
	
	
	
	def svg(self, fname = None):
		if not fname:
			fname = self.reaction.ligand
		
		# if atoms have x and y coordinates, draw them
		# otherwise draw the set of compounds
		
		# draw the thing as such
#		if self.atoms[0].x or self.atoms[0].y:
		if True:
			self.center_atoms()
			width, height = self.compute_borders()
			
			# header
			s = "<?xml version=\"1.0\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n\n<svg fill-opacity=\"1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" color-rendering=\"auto\" color-interpolation=\"auto\" stroke=\"black\" text-rendering=\"auto\" stroke-linecap=\"square\" width=\"%d\" stroke-miterlimit=\"10\" stroke-opacity=\"1\" shape-rendering=\"auto\" fill=\"black\" stroke-dasharray=\"none\" font-weight=\"normal\" stroke-width=\"1\" height=\"%d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"&apos;Dialog&apos;\" font-style=\"normal\" stroke-linejoin=\"miter\" font-size=\"12\" stroke-dashoffset=\"0\" image-rendering=\"auto\">\n  <g>\n" % (width, height)
			
			# bonds
			s += "      <g text-rendering=\"optimizeLegibility\" stroke-width=\"1.1\" shape-rendering=\"geometricPrecision\">\n"
			for b in self.bonds:
				x1 = b.source.x
				x2 = b.target.x
				y1 = b.source.y
				y2 = b.target.y
				
				s += "        <line x1=\"%f\" fill=\"none\" y1=\"%f\" x2=\"%f\" y2=\"%f\" />\n" % (x1,y1,x2,y2)
			s += "      </g>\n"
			
			# atom boxes
			s += "      <g font-size=\"15px\" fill=\"white\" text-rendering=\"optimizeLegibility\" font-family=\"Helvetica\" shape-rendering=\"geometricPrecision\" stroke=\"rgb(0,0,0)\">\n"
			for a in self.atoms:
				s += "        <rect x=\"%d\" y=\"%d\" fill=\"white\" width=\"15\" height=\"11\" stroke=\"none\" />\n" % (a.x+5,a.y-2)
#				if a.symbol == "C":
#					continue
				s += "        <rect x=\"%d\" y=\"%d\" fill=\"white\" width=\"15\" height=\"16\" stroke=\"none\" />\n" % (a.x-7,a.y-8)
			s += "      </g>\n"
			
			# atoms
			s += "      <g font-size=\"15px\" fill=\"rgb(0,0,0)\" text-rendering=\"optimizeLegibility\" font-family=\"Helvetica\" shape-rendering=\"geometricPrecision\" stroke=\"rgb(0,0,0)\">\n"
			for a in self.atoms:
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
		
		else:
			# do basic scaling
			for c in self.compoundatoms:
				self.scale_compound(c)
			# make 2-column layout
			# find the widest
			maxwidth = max([a.coords[c].x for c,atomsets in self.compoundatoms.items() for a in atomsets[0]])
			i = 0
			n = len(self.compoundatoms)
			middle = n/2
			height = 0
			width = 0
			for c in self.compoundatoms:
				if i == middle:
					width = maxwidth
					height = 0
				i += 1
				
				self.scale_compound(c, 24, width, height)
				
				print "height of", c, "is", max([a.coords[c].y for a in self.compoundatoms[c][0]])
				height = max([a.coords[c].y for a in self.compoundatoms[c][0]])
#				print "new height", height
			
			width = max([a.coords[c].x for c,atomsets in self.compoundatoms.items() for a in atomsets[0]]) + 50
			height = max([a.coords[c].y for c,atomsets in self.compoundatoms.items() for a in atomsets[0]]) + 50
			
			# header
			s = "<?xml version=\"1.0\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n\n<svg fill-opacity=\"1\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" color-rendering=\"auto\" color-interpolation=\"auto\" stroke=\"black\" text-rendering=\"auto\" stroke-linecap=\"square\" width=\"%d\" stroke-miterlimit=\"10\" stroke-opacity=\"1\" shape-rendering=\"auto\" fill=\"black\" stroke-dasharray=\"none\" font-weight=\"normal\" stroke-width=\"1\" height=\"%d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"&apos;Dialog&apos;\" font-style=\"normal\" stroke-linejoin=\"miter\" font-size=\"12\" stroke-dashoffset=\"0\" image-rendering=\"auto\">\n  <g>\n" % (width, height)
			
			for c in self.compoundatoms:
				atoms = self.compoundatoms[c][0]
				bonds = self.compoundbonds[c][0]
				
				s += "      <g>\n"
				
				# bonds
				s += "       <g text-rendering=\"optimizeLegibility\" stroke-width=\"1.1\" shape-rendering=\"geometricPrecision\">\n"
				for b in bonds:
					x1 = b.source.coords[c].x
					x2 = b.target.coords[c].x
					y1 = b.source.coords[c].y
					y2 = b.target.coords[c].y
					
					s += "         <line x1=\"%f\" fill=\"none\" y1=\"%f\" x2=\"%f\" y2=\"%f\" />\n" % (x1,y1,x2,y2)
				s += "       </g>\n"
				
				# atom boxes
				s += "       <g font-size=\"15px\" fill=\"white\" text-rendering=\"optimizeLegibility\" font-family=\"Helvetica\" shape-rendering=\"geometricPrecision\" stroke=\"rgb(0,0,0)\">\n"
				for a in atoms:
					s += "         <rect x=\"%d\" y=\"%d\" fill=\"white\" width=\"15\" height=\"11\" stroke=\"none\" />\n" % (a.coords[c].x+5,a.coords[c].y-2)
					s += "         <rect x=\"%d\" y=\"%d\" fill=\"white\" width=\"15\" height=\"16\" stroke=\"none\" />\n" % (a.coords[c].x-7,a.coords[c].y-8)
				s += "       </g>\n"
				
				# atoms
				s += "       <g font-size=\"15px\" fill=\"rgb(0,0,0)\" text-rendering=\"optimizeLegibility\" font-family=\"Helvetica\" shape-rendering=\"geometricPrecision\" stroke=\"rgb(0,0,0)\">\n"
				for a in atoms:
					s += "         <text font-size=\"11px\" fill=\"blue\" xml:space=\"preserve\" x=\"%d\" y=\"%d\" stroke=\"none\">%s</text>\n" % (a.coords[c].x+6,a.coords[c].y+8,str(a.id+1))
					
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
					
					s += "         <text fill=\"%s\" xml:space=\"preserve\" x=\"%d\" y=\"%d\" stroke=\"none\">%s</text>\n" % (color, a.coords[c].x-5,a.coords[c].y+5,a.symbol)
				s += "       </g>\n"
				s += "      </g>\n"
			
			# end
			s += "  </g>\n</svg>\n"
			
		f = open(fname, "w")
		f.write(s)
		f.close()
	
	def __str__(self):
		s = "ReactionGraph: " + str(len(self.atoms)) + " atoms, " + str(len(self.bonds)) + " bonds\n"
		
		for a in self.atoms:
			s += str(a) + "\n"
		for b in self.bonds:
			s += str(b) + str(b.reactions) + "\n"
		
		return s
	
	def __len__(self):
		return len(self.atoms)
