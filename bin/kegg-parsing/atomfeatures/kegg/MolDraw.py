#!/usr/bin/python
#
# a class to layout molecules in a pleasing manner
#



import kegg
import math, itertools, operator, sys

EPSILON = 0.00001


def pairs(l):
	p = []
	for i in range(len(l)-1):
		for j in range(i+1,len(l)):
			p.append( (l[i],l[j]) )
	return p

def equidistant_point(atoms):
	pairs = []
	for i in range(len(atoms)-1):
		for j in range(i+1,len(atoms)):
			pairs.append( (atoms[i], atoms[j]) )
	
	normals = []
	for a1,a2 in pairs:
		normals.append( normalbetween(a1,a2) )
	
	points = []
	for i in range(len(normals)-1):
		for j in range(i+1,len(normals)):
			points.append( intersection_of_normals( normals[i], normals[j] ) )
	
	avgpoint = ( sum([p[0] for p in points])/len(points) , sum([p[1] for p in points])/len(points) )
	
	return avgpoint


# intersection point between two lines
def intersection_of_normals(l1, l2):
	mp1,a1 = l1
	mp2,a2 = l2
	
	lengthx = (mp1[0] - mp2[0]) / (math.cos(a1) - math.cos(a2))
	lengthy = (mp1[1] - mp2[1]) / (math.sin(a1) - math.sin(a2))
	
	return ( mp1[0] + lengthx*math.cos(a1) , mp1[1] + lengthy*math.sin(a1) )


# normal between two atoms
def normalbetween(a1,a2):
	midpoint = ( (a1.x+a2.x)/2, (a1.y+a2.y)/2 )
	angle = getangle(a1,a2) + (math.pi/2)
	
	return midpoint, angle

# euclidean distance
def distance(a1,a2):
	try:
		x1 = a1.x
		y1 = a1.y
	except:
		x1,y1 = a1[0],a1[1]
	
	try:
		x2 = a2.x
		y2 = a2.y
	except:
		x2,y2 = a2[0],a2[1]
	
	return math.sqrt( (x1-x2)**2 + (y1-y2)**2 )

# angle between two points (zero angle = "east")
def getangle(a1, a2):
	try:
		x1 = a1.x
		y1 = a1.y
	except:
		x1,y1 = a1[0],a1[1]
	
	try:
		x2 = a2.x
		y2 = a2.y
	except:
		x2,y2 = a2[0],a2[1]
	
	a = math.atan2(y2-y1, x2-x1)  # take the angle between atoms
	while a < 0:
		a += 2*math.pi
	return a

# average angle between a set of angles
# computed by placing points on unit circle and taking their center point
# -> its angle
def avgangle(*args):
	EPSILON = 0.0000001
	
	angles = []
	for a in args:
		try:
			for elem in a:
				angles.append(elem)
		except:
			angles.append(a)
	
	# assume angles as points on unit circle, take average,
	# and solve for direction
	
	points = [(math.cos(a),math.sin(a)) for a in angles]
	avgpoint = ( sum([p[0] for p in points]) / len(points),  sum([p[1] for p in points]) / len(points) )
#	if avgpoint[0] < EPSILON and avgpoint[1] < EPSILON:
#		return angles[0] + (math.pi/2)
	
	return math.atan2(avgpoint[1], avgpoint[0])

# rotate a set of atoms by 'theta'
# optionally define a center atom to act as origo
def rotate(atoms, theta=0.0, center=None):
	if not center:
		avgx = sum([a.x for a in atoms]) / len(atoms)
		avgy = sum([a.y for a in atoms]) / len(atoms)
	else:
		avgx = center.x
		avgy = center.y
	
	for a in atoms:
		a.x, a.y = (a.x-avgx)*math.cos(theta) - (a.y-avgy)*math.sin(theta) + avgx, (a.x-avgx)*math.sin(theta) + (a.y-avgy)*math.cos(theta) + avgy


# draw ring to 'atoms' in order. The starting point is x,y, drawing angle is 'angle'
# and direction defines clockwise/counter clockwise rotation
def draw_ring(atoms):#, x=0.0, y=0.0, angle=0.0, direction=1):
	# assigned atoms, don't change them
	assigned = len( filter(lambda x:x is not None, [a.x for a in atoms]) )
	k = len(atoms)
	
	if assigned == 0:
		coords = ringcoords(k)
		for i in range(k):
			atoms[i].x = coords[i][0]
			atoms[i].x = coords[i][1]
	# one of the atoms is already fixed, use that as the starting point
	elif assigned == 1:
		# find the fixed atom
		a = [a for a in atoms if a.x is not None][0]
		fixedneighs = [ne for ne in a.GetAtomNeighbors() if a.x is not None]
		angles_to_fixeds = [getangle(a,ne) for ne in fixedneighs]
		outangle = avgangle(angles_to_fixeds)
		inangle = outangle + math.pi
		# sort to start from fixed atom
		while atoms[0] != a:
			atoms = atoms[1:] + [atoms[0]]
		coords = ringcoords(k, x=a.x, y=a.y, angle=inangle + 2*math.pi/len(atom)/2, direction=-1)
		for i in range(k):
			atoms[i].x = coords[i][0]
			atoms[i].x = coords[i][1]
	# two of the atoms are already fixed, use those as starting points
	elif assigned == 2:
		fixed_atoms = [a for a in atoms if a.x is not None]
		ca = fixed_atoms[0] # center (shared) atom
		sa = fixed_atoms[1] # other shared atom
		# next atom from sa in fixedatoms
		current_ring = list(reduce(set.intersection, [set(self.supernodes[a]) for a in atoms]))[0]
		outer_ring = list(set(self.supernodes[sa]) & set(self.supernodes[ca]) - set([current_ring]))[0]
		next_outer = [x for x in sa.GetAtomNeighbors() if outer_ring in self.supernodes[x] and x != ca][0]
		
		outangle = getangle(ca,sa)
		nextangle = getangle(sa, next_outer)
		if 0 < nextangle - outangle < math.pi:
			direction = -1
		else:
			direction = 1
		
		# put centeratom in first place
		while ca != atoms[0]:
			atoms = atoms[1:] + [atoms[0]]
		# put other shared atom at second place
		if sa != atoms[1]:
			atoms = [atoms[0]] + list(reversed(atoms[1:]))
		
		coords = ringcoords(k, 25, ca.x, ca.y, outangle, direction)
		for i in range(k):
			atoms[i].x = coords[i][0]
			atoms[i].y = coords[i][1]
	# three atoms are already fixed, they have to belong to different rings
	# -> we need to assign the rest 'k-3' atoms evenly
	# if we have A -> B -> C (already fixed),
	# then lets start from C and go on for 2 more atoms in case of k=5
	elif assigned == 3:
		pass
		
		
	


def ringcoords(k, bondlen=25, x=0.0, y=0.0, angle=0.0, direction=1):
	results = []
	currentx = x
	currenty = y
	for i in range(k):
		results.append( (currentx, currenty) )
		currentx += math.cos(angle + direction * i * 2 * math.pi / k) * bondlen
		currenty += math.sin(angle + direction * i * 2 * math.pi / k) * bondlen
	return results













class MolDraw:
	def __init__(self, mol):
		self.molecule = mol
		self.atoms = self.molecule.graph.atoms
		self.bonds = self.molecule.graph.bonds
		self.supergraph = kegg.Graph()  # superatoms
		self.supernodes = {}  # atoms -> superatoms
		
		self.bondlen = 20.0
		
		
	
	
	def cyclic(self):
		# test whether there exists rings of sizes [4...10]
		# start to walk from each node and check whether k'th node in sequence is starting node
		
		rings = self.get_rings()
		
		if rings:
			return True
		return False
	
	def acyclic(self):
		return not self.cyclic()
	
	
	def get_rings(self):
		ringsets = set()
		def dfs(seq, k):
			if k > 7:
				return
			
			for na in seq[-1].GetAtomNeighbors():
				if na == seq[0]:
					if 4 <= k:
						ringsets.add(tuple(sorted(seq)))
				if na not in seq:
					dfs(seq + [na], k+1)
		
		for a in self.atoms:
			dfs([a],1)
		
		# have all unique rings, restore order to rings
		rings = []
		for r in ringsets:
			rings.append([])
			start = list(r)[0]
			rings[-1].append(start)
			next = start
			
			while True:
				neighs = [a for a in next.GetAtomNeighbors() if a in r and a not in rings[-1]]
				if neighs:
					next = neighs[0]
					rings[-1].append(next)
				else:
					break
		
		# remove supersets
		while True:
			for r1,r2 in [(r1,r2) for r1 in rings for r2 in rings if r1 != r2]:
				if len(set(r1) & set(r2)) == len(r1):
					rings.remove(r2)
					break
				elif len(set(r1) & set(r2)) == len(r2):
					rings.remove(r1)
					break
			else:
				break
		
		
		
#		print "found", len(rings), "rings"
#		for r in rings:
#			print map(str, r)
		
		return rings
	
	
	
	
	
	
	
	# main layout class
	def layout(self):
		
#		for a in self.atoms:
#			print a, map(str, a.neighs)
		
		if self.cyclic():
			print "cyclic molecule"
			self.cyclic_layout()
		else:
			print "acyclic molecule"
			self.acyclic_layout()
		
		self.normalize_coordinates()
	
	
	# make top-left corner (0,0) point
	# ideally internal coordinates now correspond to external coordinates
	def normalize_coordinates(self):
		left = top = 50
		# make all coordinates positive and left-top corner -most atom into (50,50) coordinate
		min_x = min([a.x for a in self.atoms])
		max_x = max([a.x for a in self.atoms])
		min_y = min([a.y for a in self.atoms])
		max_y = max([a.y for a in self.atoms])
		
		for a in self.atoms:
			a.x -= (min_x - left)
			a.y -= (min_y - top)
	
	
	def initialize(self):
		kegg.Node.idcounter = 1
		# construct the a graph for superatoms
		for a in self.atoms:
#			newnode = self.supergraph.addnode()
#			newnode.atoms = [a]
#			newnode.ring = False
#			newnode.chain = False
#			newnode.junction = False
#			self.supernodes[a] = newnode
			
			a.x = None
			a.y = None
#		for b in self.bonds:
#			self.supergraph.addedgeto(self.supernodes[b.source], self.supernodes[b.target])
		
		self.supergraph.rings = []
		self.supergraph.junctions = []
		self.supergraph.chains = []
	
	
	def cyclic_layout(self):
		# the layout goes in phases and is based on groups with internal coordinates
		# 
		# (1) find all rings (also fused) and give them internal coordinates, and contract them into groups
		# (2) start form rings: the next atom can be put into several angles from ring (prefer most distant angle)
		#     continue and give options for each node
		#     finally optimize such that the combination of angles maximizes intra-distances
		#     
		# (last) rotate the whole molecule such that x-axis is maximized (width) 
		
		# we have:
		# - self.molecule and self.atoms for the usual atom graph
		# - supergraph which contains superatoms containing [1,...] atoms
		# - supernodes containing mapping from atoms to supernodes
		
		# initialize supergraph
		self.initialize()
		
		# mark the areas of the molecule belonging to rings, chains and junctions
		self.mark_rings()
		self.mark_chains()
		self.mark_junctions()
		
		# connect the supergraph
		self.connect_supergraph()
		
		# give internal fixed coordinates to chains and rings
		self.coordinate_rings()
		
		sys.exit(1)
		
		self.coordinate_junctions()
		self.coordinate_chains()
		
		
		# contract overlapping rings into a single ring system
		self.contract_fuserings()
		
		# set atoms to point to single supernodes now
		for a in self.atoms:
			self.supernodes[a] = self.supernodes[a][0]
		
	
		# start the layout algorithm:
		# 1) start from initial ring, fix its location
		# 2) go through the neighbors of the initial ring:
		#    set their angles to "optimal" angles based on degree
		#    d=3 -> 120 degrees
		#    d=4 -> 90 degrees
		# 3) go through the neighbors of the processed region, do the same
		#
		
#		for a in self.atoms:
#			print a, a.x, a.y, self.supernodes[a]
		
		def degreesort( (c1,n1), (c2,n2) ):
			if len(n1.neighs) > len(n2.neighs):
				return 1
			elif len(n1.neighs) < len(n2.neighs):
				return -1
			return 0
			
		
		print
		print "before layouting internal coords are..."
		for a in self.atoms:
			print "%d: (%d,%d)" % (a.id, int(a.x), int(a.y))
		
		current = sorted([x for x in self.supergraph.nodes if x.ring], key=lambda x: len(x.atoms))[-1]
		closed = [current]
		closedatoms = current.atoms
		border = [(current,ne) for ne in current.neighbors()]
		border.sort(degreesort)
		
		print
		print "starting layout algorithm from largest ringsystem"
		print current, map(str, current.atoms)
		for a in current.atoms:
			print "", a.id, a.symbol, a.x, a.y
		
		
		while border:
			current,next = border.pop()
			
			if next in closed:
				continue
			
			centeratom,nextatom = [(a1,a2) for a1 in current.atoms for a2 in next.atoms if a1.IsNeighbor(a2)][0]
			neighs = [x for x in centeratom.GetAtomNeighbors() if x in closedatoms]
			outneighs = [x for x in centeratom.GetAtomNeighbors() if x not in closedatoms]
			
#			if centeratom.Degree() == 2:
#				outangle = math.pi + getangle(centeratom, neighs[0])
#			elif centeratom.Degree() == 3:
#				outangle = 2*math.pi/3 + getangle(centeratom, ne)
#			elif centeratom.Degree() == 4:
#				if len(neighs) == 1:
#					outangle = math.pi + getangle(centeratom, neighs[0])
#				elif len(neighs) == 2:
#					outangle = math.pi/2 + getangle(centeratom, neighs[0])
#				elif len(neighs) == 3:
#					outangle = math.pi + avgangle(angles)
			
			
			angles = []
			for ne in neighs:
				angles.append( getangle(centeratom, ne) )
				print "neighbor", ne, "with angle", angles[-1]
			
			innerdir = avgangle(angles)
			outerdir = math.pi + innerdir
			
#			if len(outneighs) == 2:
#				
			
#			outerdir = centeratom.freeangles.pop(0) + getangle(centeratom, neighs[0])
#			nextatom.freeangles.pop(0)
			
#			if len(neighs) == 1:
#				outerdir = int(centeratom.Degree()/2) * 2*math.pi/centeratom.Degree() + getangle(centeratom, neighs[0])
#			if len(neighs) == 2:
#				outerdir = centeratom.freeangles[0]
#				del centeratom.freeangles[0]
#				
#				outdirs = [2*math.pi/centeratom.Degree() + getangle(centeratom, ne) for ne in neighs]
#				
#				outerdir = 2*math.pi/centeratom.Degree() + getangle(centeratom, neighs[0])
			
			x = centeratom.x + math.cos(outerdir) * self.bondlen
			y = centeratom.y + math.sin(outerdir) * self.bondlen
			
			print "centeratom:", centeratom.x, centeratom.y, "outangle:", outerdir, "innerangle:", innerdir
			self.calibrate(next, nextatom, centeratom, x, y)
			
			for ne in next.neighbors():
				if ne not in border and ne not in closed:
					border.append((next,ne))
			
			closed.append(next)
			for a in next.atoms:
				closedatoms.append(a)
			border.sort(degreesort)
		
		print
		print "after layout coords are..."
		for a in self.atoms:
			print a, a.x, a.y
		
		self.rotate_wide()
		self.normalize_coordinates()
#		for a in self.atoms:
#			print a, a.x, a.y
		
	
	def rotate_wide(self):
#		print "widening..."
		# rotate so that the structure is as wide as possible
		oldwidth = 0
		newwidth = max([a.x for a in self.atoms]) - min([a.x for a in self.atoms])
		
		while newwidth > oldwidth:
#			print oldwidth, newwidth
			oldwidth = newwidth
			rotate(self.atoms, 0.01)
			newwidth = max([a.x for a in self.atoms]) - min([a.x for a in self.atoms])
		
		oldwidth = 0
		newwidth = max([a.x for a in self.atoms]) - min([a.x for a in self.atoms])
		
		while newwidth > oldwidth:
#			print oldwidth, newwidth
			oldwidth = newwidth
			rotate(self.atoms, -0.01)
			newwidth = max([a.x for a in self.atoms]) - min([a.x for a in self.atoms])
		
	
	
	def calibrate(self, node, atom, prevatom, x, y):
		# we calibrate the 'node' region, starting from 'atom', going from 'prevatom'
		#'atom' is given a new coordinate, calibrate rest of the node
		
		# translate
		dx = x - atom.x
		dy = y - atom.y
		
		print "calibrating", map(str, node.atoms), "with", dx, dy, x, y
		
		atom.x = x
		atom.y = y
		
		print "", atom, "changed to", x,y
		for a in node.atoms:
			if a != atom:
				a.x += dx
				a.y += dy
				print "", a, "changed to", a.x, a.y
		
		# get angle from prevatom -> atom
		#           from atom -> centeraotm
		# to be the same
		
		centeratom = kegg.Atom("X")
		centeratom.x = sum([a.x for a in node.atoms]) / len(node.atoms)
		centeratom.y = sum([a.y for a in node.atoms]) / len(node.atoms)
		
		angle1 = getangle(prevatom, atom)
		angle2 = getangle(atom, centeratom)
		anglediff = angle1 - angle2
		
		# rotate
		rotate(node.atoms, anglediff, atom)
		
		print "after rotation:"
		for a in node.atoms:
			print "", a, a.x, a.y
		
	
	
	
	def coordinate_chains(self):
		print "fixing chain's internal coordinates..."
		
		for n in self.supergraph.nodes:
			if n.chain:
				if len(n.atoms) == 1:
					n.atoms[0].x = 0.0
					n.atoms[0].y = 0.0
				elif len(n.atoms) == 2:  # 2-long chain -> put in straight line
					n.atoms[0].x, n.atoms[0].y = 0.0, 0.0
					n.atoms[1].x, n.atoms[1].y = self.bondlen, 0.0
				else:  # longer chain -> zigzag pattern
					sx, sy = 0.0, 0.0
					angle = math.pi/6  # 30 degrees
					for a in n.atoms:
						a.x = sx
						a.y = sy
						sx += math.cos(angle) * self.bondlen
						sy += math.sin(angle) * self.bondlen
						angle *= -1 # flip angle
				print " chain [%s] coords: [%s]" % ( ",".join([str(a.id) for a in n.atoms]), ",".join([ "("+str(int(a.x))+","+str(int(a.y))+")" for a in n.atoms ])  )
			
	
	def coordinate_junctions(self):
		print "fixing junctions's internal coordinates"
		for n in self.supergraph.nodes:
			if n.junction:
				n.atoms[0].x = 0.0
				n.atoms[0].y = 0.0
			print " junction %d coords (0,0)" % (n.atoms[0].id)
	
	
	
	def coordinate_rings(self):
		print "fixing ring's internal coordinates..."
		
		# largest rings first
		rings = [n for n in self.supergraph.nodes if n.ring]
		rings.sort(key=lambda x: len(x.atoms), reverse=True)
		i = 1
		for r in rings:
			r.ringid = i
			i += 1
		
		# find ringsystems (connected rings)
		ringsystems = []
		border = []
		closed = []
		ringmap = {}
		for r in rings:
			if r in closed:
				continue
			
			ringsystems.append( kegg.Graph() )
			ringgraph = ringsystems[-1]
			rn = ringgraph.addnode()
			rn.ring = r
			closed.append(r)
			ringmap[r] = rn
			
			for ne in r.neighbors():
				if ne.ring:
					border.append( (r,ne) )
			
			while border:
				prev,next = border.pop()
				
				if next in ringmap:
					ringgraph.addedgeto(ringmap[prev], ringmap[next] )
				else:
					rn = ringgraph.addnodeto(ringmap[prev])
					rn.ring = next
					closed.append(next)
					ringmap[next] = rn
				
				for ne in next.neighbors():
					if ne.ring and ne not in closed:
						border.append( (next,ne) )
		
		
		
		self.seqs = []
		def dfs(sequence, border, n):
			print "dfs", map(lambda x:x.ringid, sequence), map(lambda x:x.ringid,border), n.ringid,
#			print map(str, sequence), map(str, border), n
			
			if not sequence:  # first element
				self.coordinate_ring(n.atoms)
				success = True
			else:
				success = self.coordinate_fuse_ring(n)
			
			if not success:
				print "fail" 
				# save result
				if not self.seqs or len(sequence) == self.seqs[0]:
					self.seqs.append(sequence)
				elif len(sequence) > self.seqs[0]:
					self.seqs = [sequence]
				return
			
			print
			
			sequence.append(n)  # update sequence
#			print "adding", n
			
			for ne in n.neighbors():  # update border
				if ne.ring and ne not in sequence and ne not in border:
					border.append(ne)
			
			# border empty
			if not border:
				# save result
				if not self.seqs or len(sequence) == self.seqs[0]:
					self.seqs.append(sequence)
				elif len(sequence) > self.seqs[0]:
					self.seqs = [sequence]
			else:
				for next in border:
					dfs(sequence, border, next)
		
		
		for n in self.supergraph.nodes:
			print n, map(str, n.neighs), map(str, n.atoms)
		
		
		for r in rings:
			print r.ringid, map(str, r.atoms)
		
		print "starting dfs"
		for n in rings:
			dfs([], [], n)
		
		for s in self.seqs:
			print map(str,s)
		sys.exit(1)
		
		
		closed = []
		for n in rings:
			if n in closed:
				continue
			
			self.coordinate_ring(n.atoms)
			closed.append(n)
			print " ring %d:[%s] coords [%s]" % ( n.ringid, ",".join([str(a.id) for a in n.atoms]), " ".join(["(%d,%d)" % (a.x,a.y) for a in n.atoms])  )
			
			# all neighboring rings to border
			border = [neigh for neigh in n.neighbors() if neigh not in closed and neigh.ring]
			while border:
				next = border.pop(0)
				if next in closed:
					continue
				
				fixedatoms = [a for a in next.atoms if a.x is not None]
				
				if len(fixedatoms) == 0:
					continue  # don't handle if neighboring ring is attached by covalent bond
				else:
					result = self.coordinate_fuse_ring(next)
					# the ring didn't fit -> remove that ring
					if result == False:
						next.ring = False
						print "didn't work"
					else:
						print " -> ring %d:[%s] coords [%s]" % (next.ringid, ",".join([str(a.id) for a in next.atoms]), " ".join(["(%d,%d)" % (a.x,a.y) for a in next.atoms])  )
					closed.append
				
				for ne in next.neighbors():
					if ne not in closed and ne not in border and ne.ring:
						border.append(ne)
				
	
	# coordinate a ringnode when at least one atom is already fixed
	def coordinate_fuse_ring(self, ringnode):
		atoms = ringnode.atoms
		fixedatoms = [a for a in atoms if a.x is not None]
		k = len(atoms)
		
		# a single shared atom, direct this ring outwards from the shared-atom
		if len(fixedatoms) == 1:
			# lets check shared-atoms neighbors
			centeratom = fixedatoms[0]
			outerneighs = [a for a in centeratom.GetAtomNeighbors() if a.x is not None]
			angles = [getangle(centeratom,a) for a in outerneighs]
			
			outangle = avgangle(angles)
			inangle = math.pi + outangle
			
			while atoms[0] != centeratom: # rotate until fixed atom is first
				atoms = atoms[1:] + [atoms[0]]
			
			self.coordinate_ring(atoms, inangle - 2*math.pi/k)
		
		# two shared atoms -> a shared edge
		elif len(fixedatoms) == 2:
			ca = fixedatoms[0] # center (shared) atom
			sa = fixedatoms[1] # other shared atom
			next = [a for a in sa.GetAtomNeighbors() if a.x is not None and a != ca][0]  # next atom from sa outside
			
			outangle  = getangle(ca, sa)
			nextangle = getangle(sa, next)
			if 0 < nextangle - outangle < math.pi:
				direction = -1
			else:
				direction = 1
			
			# put centeratom in first place
			while ca != atoms[0]:
				atoms = atoms[1:] + [atoms[0]]
			# put other shared atom at second place
			if sa != atoms[1]:
				atoms = [atoms[0]] + list(reversed(atoms[1:]))
			
			self.coordinate_ring(atoms, outangle, ca.x, ca.y, direction)
		
		# if 3 consecutive fixed atoms -> 
		elif len(fixedatoms) == 3:
			consecutive = False
			for i in range(k):
				if atoms[i].x is not None:
					if atoms[i+1 % k].x is not None and atoms[i+2 % k].x is not None:
						consecutive = True
						break
			
			if not consecutive:
				print "error, fixed atoms are not consecutive"
				return False
			
			# wind fixed atoms to first place
			while atoms[0].x is None:
				atoms = atoms[1:] + [atoms[0]]
			
			# get middle point
			centerpoint = equidistant_point( atoms[0:3] )
			radius = distance(centerpoint, (atoms[0].x,atoms[0].y))
			
			alpha = 2*math.pi - (getangle(centerpoint, atoms[0]) - getangle(centerpoint, atoms[2]))
			angle = alpha / (k-2)
			centerpoint + ( radius*math.cos(angle), radius*math.sin(angle) )
			
			startangle = getangle(centerpoint, atoms[2])
			if getangle(centerpoint, atoms[2]) > getangle(centerpoint, atoms[1]):
				angle *= -1
			
			currangle = startangle + angle
			
			for a in atoms[3:]:
				a.x = centerpoint[0] + radius*math.cos(currangle)
				a.y = centerpoint[1] + radius*math.sin(currangle)
				currangle += angle
		
		else:
			# remove ring status
			return False
		
		return True
	
	
	def contract_fuserings(self):
		print "contracting fuserings"
		ringnodes = [n for n in self.supergraph.nodes if n.ring]
#		print map(str,ringnodes)
		while True:
			for r1,r2 in [(r1,r2) for r1 in ringnodes for r2 in ringnodes if r1 != r2]:
#				print r1,r2, set(r1.atoms) & set(r2.atoms)
				if set(r1.atoms) & set(r2.atoms):
#					print r1,r2
					r1.atoms += r2.atoms
					r1.atoms = list(set(r1.atoms))
					
					for a in r2.atoms:
						self.supernodes[a] = [r1]
					
					for ne in self.supergraph.neighbors([r2]):
						if ne not in r1.neighs:
							self.supergraph.addedgeto(r1,ne)
					
					self.supergraph.removenode(r2)
					ringnodes.remove(r2)
#					print "contracted ring", r2
					
					break
			else:
				break
	
	
	
	# place coordinates for rings
	# start form 'startx,starty'
	# next one is towards direction 'angle'
	# after that counter/clockwise depending on 'direction'
	def coordinate_ring(self, atoms, angle=0.0, startx=0.0, starty=0.0, direction=1):
		k = len(atoms)
		sx = startx
		sy = starty
		for i in range(k):
			a = atoms[i]
			a.x = sx
			a.y = sy
			
			# dir defines counter/clockwise path
			sx += math.cos(angle + direction * (i) * 2 * math.pi / k) * self.bondlen
			sy += math.sin(angle + direction * (i) * 2 * math.pi / k) * self.bondlen
	
	
	
	def mark_chains(self):
		# mark sequences of atoms...
		# chain is a connected sequence of nodes of degree 2 or 1
		
		print "finding chains..."
		
		chains = []
		closed = []
		# start from each atom
		for n in self.atoms:
			# if already processed, or ring or junction, or degree too large -> pass
			if n in closed or n in [a for a in self.supergraph.rings] or n in [a for a in self.supergraph.junctions] or n.Degree() not in (1,2):
				continue
			
			chain = [n]
			neighs = n.GetAtomNeighbors()
			closed.append(n)
			
			# iteratively go through the border of the chain, adding to if degree = 1 or 2
			while neighs:
				curr = neighs.pop()
				
				if curr.Degree() in (1,2) and curr not in closed+neighs and curr not in [a for a in self.supergraph.rings] and curr not in [a for a in self.supergraph.junctions]:
					chain.append(curr)
					closed.append(curr)
					
					for ne in curr.GetAtomNeighbors():
						neighs.append(ne)
			
			chains.append(chain)
		
		# create supernodes for chains
		for chainatoms in chains:
#			chainnodes = []  # list(set(reduce(operator.add, [self.supernodes[a] for a in chainatoms])))
#			chainborder = self.supergraph.neighbors( chainnodes )
			
			chainnode = self.supergraph.addnode()
			chainnode.atoms = chainatoms
			chainnode.chain = True
			chainnode.ring = False
			chainnode.junction = False
			self.supergraph.chains.append(chainnode)
			
#			for bn in chainborder:
#				self.supergraph.addedgeto(chainnode, bn)
#			for n in chainnodes:
#				self.supergraph.removenode(n)
#			for cn in chainatoms:
#				self.supernodes[cn] = [chainnode]
		
#		self.supergraph.chains = chains
		
		print " found %d chains" % (len(chains))
		
		
#	def connect_supergraph(self):
#		for j,c in pairs(self.supergraph.junctions, self.supergraph.chains):
#			if ...
#			pass
	
	
	
	
	def mark_junctions(self):
		print "finding junctions..."
		
		juncs = []
		for a in self.atoms:
			if a.Degree() > 2 and a not in [a for r in self.supergraph.rings for a in r.atoms]:
#				self.supernodes[a][0].junction = True
				newnode = self.supergraph.addnode()
				newnode.atoms = [a]
				newnode.chain = False
				newnode.ring = False
				newnode.junction = True
				self.supergraph.junctions.append(newnode)
				juncs.append(newnode)
		
		for j1,j2 in pairs(juncs):
			if j1.atoms[0].IsNeighbor(j2.atoms[0]):
				self.supergraph.addedgeto(j1,j2)
		
		print " found %d junctions" % (len(juncs))
	
	
	
	# add supernodes in middle of rings
	def mark_rings(self):
		
		print "finding rings..."
		
		rings = self.get_rings()
		
#		removables = []
		ringnodes = []
		for ring in rings:
#			print map(str, ring)
			# add new superatom as a component
			newnode = self.supergraph.addnode()
			newnode.atoms = ring
			newnode.ring = True
			newnode.chain = False
			newnode.junction = False
			ringnodes.append(newnode)
		
		# connect rings
		for rn1,rn2 in pairs(ringnodes):
			if set(rn1.atoms) & set(rn2.atoms):
				self.supergraph.addedgeto(rn1, rn2)
			else:
				if set(rn1.atoms) & set([ne for a in rn2.atoms for ne in a.GetAtomNeighbors()]):
					self.supergraph.addedgeto(rn1, rn2)
		
		self.supergraph.rings = ringnodes
		
#			neighs = set([ne for a in ring for ne in self.supernodes[a].neighbors() if ne != newnode])
#			for n in neighs:
#				self.supergraph.addedgeto(newnode, n)
#			
#			for a in ring:
#				if not self.supernodes[a].ring:
#					removables.append(self.supernodes[a])
#				self.supernodes[a] = newnode
		
		
#		removables = set(removables)
		
#		print map(str, removables)
		
		# connect ringnodes
#		for b in self.bonds:
#			end1 = self.supernodes[b.source]
#			end2 = self.supernodes[b.source]
#			if end1 != end2 and end1.ring and end2.ring:
#				self.supergraph.addedgeto(end1,end2)
		
#		for n in self.supergraph.nodes:
#			print n, map(str, n.atoms)
		
		# remove ring's nodes
#		for rem in removables:
#			print "removing", rem
#			self.supergraph.removenode(rem)
		
#		for n in self.supergraph.nodes:
#			print n, map(str, n.atoms)
		
#		for a in self.atoms:
#			self.supernodes[a] = []
#		for rn in ringnodes:
#			for a in rn.atoms:
#				self.supernodes[a].append(rn)
#		for n in self.supergraph.nodes:
#			if len(n.atoms) == 1:
#				self.supernodes[n.atoms[0]] = [n]
		
		print " found %d rings" % len(rings)
		
		
		
	
	
	
	
	def acyclic_layout(self):
		pass
	

































