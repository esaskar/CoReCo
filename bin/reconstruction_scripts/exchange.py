#!/usr/bin/env python

import sys
import libsbml



class Reaction:
	def __init__(self, id, name, lb, ub, reactants, products):
		self.id = id
		self.name = name
		self.lb = lb
		self.ub = ub
		self.reactants = set(reactants)
		self.products = set(products)

	def get_id(self):
		return self.id

	def get_name(self):
		return self.name

	def get_lb(self):
		return self.lb

	def get_ub(self):
		return self.ub

	def get_reactants(self):
		return self.reactants
	
	def get_products(self):
		return self.products

class Compound:
	def __init__(self, id, name, stoichiometry):
		self.id = id
		self.name = name
		self.stoichiometry = stoichiometry

	def get_id(self):
		return self.id
	
	def get_name(self):
		return self.name
	
	def get_stoichiometry(self):
		return self.stoichiometry


class Exchange:	
	def read_reactions(self,reactions_path):
		reaction_id = ""
		reaction_name = ""
		lb = ""
		ub = ""
		is_reactant = False;
		is_product = False;	
		reactants = []
		products = []
		reactions = []

		f = open(reactions_path)
		for s in f:		
			if "Ex_" in s:				
				reaction_id, reaction_name, lb, ub = s.split(" - ");
			elif "Reactants:" in s:
				is_reactant = True
			elif is_reactant:
				if "Products:" in s:
					is_reactant = False
				else:
					if "[" in s:
						stoichiometry, reactantId, reactantName = s.split(" - ");
						reactant = Compound(reactantId, reactantName, stoichiometry)
						reactants.append(reactant)
			elif "Products:" in s:	
				is_reactant = False				
				is_product = True
			elif is_product:
				if "------------------" in s:
					reaction = Reaction(reaction_id, reaction_name, lb, ub, reactants, products)
					reactions.append(reaction)
					products = []
					reactants = []	
					is_product = False
				else:	
					if "[" in s:
						stoichiometry, productId, productName = s.split(" - ");
						product = Compound(productId, productName, stoichiometry)
						products.append(product)
			elif "-------------" in s:
				reaction = Reaction(reaction_id, reaction_name, lb, ub, reactants, products)
				reactions.append(reaction)
				products = []
				reactants = []
		f.close()	
		return reactions
				
		
	def contains_specie(self, reaction, model):			
		for reactant in reaction.get_reactants():
			isHere = False
			for specie in model.getListOfSpecies():
				if specie.getId() in reactant.get_id():
					isHere = True
					break
			if not isHere:
				return False	
		for product in reaction.get_products():
			isHere = False
			for specie in model.getListOfSpecies():
				if specie.getId() in product.get_id():
					isHere = True	
					break
			if not isHere:
				return False	
					
		return True
			
	def add_exchange_reactions(self, model_path):
		DEFAULT_COMPARTMENT = "cytosol"
		reader = libsbml.SBMLReader()
		document = reader.readSBML(model_path)
		model = document.getModel()	
		print model.getLevel()
		for reaction in self.reactions:
			## if model contains reactants and products
			if self.contains_specie(reaction, model):
				r = model.createReaction()
				r.setId(reaction.get_id())
				r.setMetaId("meta_%s" % (reaction.get_id()))	
				r.setName(reaction.get_name())
				if float(reaction.get_lb()) < 0 and float(reaction.get_ub()) > 0:
					r.setReversible(True)
				else:
					r.setReversible(False)
				r.setSBOTerm("SBO:0000176")
				r.setFast(False)
				r.setCompartment(DEFAULT_COMPARTMENT)
				law = r.createKineticLaw()
				if model.getLevel()==2:
					lbParameter = law.createParameter()
					ubParameter = law.createParameter()
				else:
					lbParameter = law.createLocalParameter()
					ubParameter = law.createLocalParameter()
				lbParameter.setId("LOWER_BOUND")
				lbParameter.setValue(float(reaction.get_lb()))
				ubParameter.setId("UPPER_BOUND")
				ubParameter.setValue(float(reaction.get_ub()))

				for reactant in reaction.get_reactants():
					re = r.createReactant()
					re.setSpecies(reactant.get_id())
					re.setStoichiometry(float(reactant.get_stoichiometry()))
					re.setConstant(True)

				for product in reaction.get_products():
					pr = r.createProduct()
					pr.setSpecies(product.get_id())
					pr.setStoichiometry(float(product.get_stoichiometry()))
					pr.setConstant(True)		
		libsbml.writeSBMLToFile(document, model_path)
		
	def __init__(self, model_path, reactions_path):
		self.reactions = []		
		self.model_path = model_path
		self.reactions = self.read_reactions(reactions_path)
		self.add_exchange_reactions(model_path)	

class Bounds:

	def setBounds(self, bounds_path, model_path):
		reader = libsbml.SBMLReader()
		document = reader.readSBML(model_path)
		model = document.getModel()		 
		f = open(bounds_path)   
		for s in f: 
			name, idreaction, model_name, lb, ub = s.strip().split("\t") 	
			name = name.replace("-", "")
		   	reaction = model.getReaction(name)
			if reaction is None:
				# test if there is a reaction starting with a
				# number, N was added to these in network2sbml.py
				reaction = model.getReaction("N"+name)
		    	if reaction is not None:
				law = reaction.createKineticLaw()
				if model.getLevel()==2:
					lbParameter = law.createParameter()
					ubParameter = law.createParameter()
				else:
					lbParameter = law.createLocalParameter()
					ubParameter = law.createLocalParameter()
				lbParameter.setId("LOWER_BOUND")
				lbParameter.setValue(float(lb))
				ubParameter.setId("UPPER_BOUND")
				ubParameter.setValue(float(ub))
				if(float(ub)>0 and float(lb)<0):
					reaction.setReversible(True)
				else:
					reaction.setReversible(False)
		f.close()		
		libsbml.writeSBMLToFile(document, model_path)

	def __init__(self, bounds_path, model_path):
		self.setBounds(bounds_path, model_path)
class Formulas:

	def setFormulas(self, formulas_path, model_path):
		reader = libsbml.SBMLReader()
		document = reader.readSBML(model_path)
		model = document.getModel()		 
		f = open(formulas_path)   
		for s in f: 
			if not s.startswith("#"):
				molid, formula = s.strip().split("\t") 	
				species = model.getSpecies(molid)
				if species is not None:
					note = """
<notes>
  <body xmlns="http://www.w3.org/1999/xhtml">
    <p>FORMULA:%s</p>
  </body>
</notes>
""" % (formula.replace(":","").replace(",",""))
					if species.isSetNotes():
						species.appendNotes(note)
					else: 
						species.setNotes(note)
				
					if not species.isSetNotes():
						print "Failed to assign note to %s:\n%s" % (molid,note)

		
		f.close()		
		libsbml.writeSBMLToFile(document, model_path)

	def __init__(self, formulas_path, model_path):
		self.setFormulas(formulas_path, model_path)

if __name__ == "__main__":
 	import sys
	Exchange(sys.argv[1], sys.argv[2])
	Bounds(sys.argv[3], sys.argv[1])
	if len(sys.argv) > 4:
		print "Setting Formulas"
		Formulas(sys.argv[4], sys.argv[1])
		
 
	
