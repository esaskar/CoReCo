#!/usr/bin/python
# -*- coding: iso-8859-15 -*-
#
# Message-classes
# Message types are real value, boolean value, set, distribution and vector
#
# currently only vector-type is used, other styles could also be used
# for example the ring detection could work by using as msgs:s distance vectors and
# as node values booleans. Then the boolean would be true if the distance vector has 6 
# in the correct position and msg vectors are accumulated through iterations.
#
#


class Real:
	def __init__(self, value=0.0, node=None):
		self.value = value
		self.node = node

	def update(self, other):
		self.value += other.value
	
	def __div__(self, other):
		return Real(self.value / other.value, self.node)
	
	def __add__(self, other):
		if isinstance(other, Real):
			return Real(self.value + other.value, self.node)
		else:
			return Real(self.value + other, self.node)
	
	def __str__(self):
		return str(self.value)
	

class Set:
	pass

class Distribution:
	pass

class Boolean:
	pass

class Vector:
	# represents a vector of type [0 0 ... val ... 0 0], where
	# for each node its value is nonzero and others are zero
	
	def __init__(self, size, default=None, position=0, value=None, node=None):
		self.vec = [default for x in range(size)]
		self.vec[position] = value
		self.size = size
		self.node = node
	
	def fill(self, value):
		self.vec = [value for x in range(self.size)]
		return self
	
	def update(self, other):
		self.vec = [ min(self.vec[i], other.vec[i]) for i in range(self.size) ]
#		for i in range(0, self.size):
#			self.vec[i] = self.vec[i] + other.vec[i]
	
	def sum(self):
		vals = [v for v in self.vec if v is not None]
		return sum(vals)
	
	# sets all values to None which are larger than in 'other'
	def min_or_none(self, other):
		self.vec = [ (self.vec[i] if other.vec[i] is None or self.vec[i] <= other.vec[i] else None) for i in range(self.size) ]
		
	
	def min(self, other):
		for ind,val in enumerate(self.vec):
			if self.vec[ind] is None:
				self.vec[ind] = other.vec[ind]
			elif other.vec[ind] is None:
				self.vec[ind] = self.vec[ind]
			else:
				self.vec[ind] = min(self.vec[ind], other.vec[ind])
		
#		self.vec = [ min(self.vec[i], other.vec[i]) for i in range(self.size) ]
	
	def max(self, other):
#		print "z", self.vec, other.vec
		self.vec = [ max(self.vec[i], other.vec[i]) for i in range(self.size) ]
#		print self.vec
	
	def __add__(self, other):
		
#		print "Z", type(other), other
#		
		if type(other) == int:
			other = Vector(self.size).fill(1)
		
		v = Vector(self.size)
		v.node = self.node
		for i in range(0,self.size):
			if self.vec[i] is None or other.vec[i] is None:
				v.vec[i] = None
			else:
				v.vec[i] = self.vec[i] + other.vec[i]
		return v
	
	def __sub__(self, other):
		v = Vector(self.size)
		v.node = self.node
		for i in range(0,self.size):
			v.vec[i] = self.vec[i] - other.vec[i]
		return v
	
	def __div__(self, other):
		v = Vector(self.size)
		for i in range(0,self.size):
			v.vec[i] = self.vec[i] / other.vec[i]
		return v
	
	def __eq__(self, other):
		return self.vec == other.vec
	
	def __ne__(self, other):
		return self.vec != other.vec
	
	def __str__(self):
		return str(self.vec)

