#!/usr/bin/env python


import libsbml


class Reaction:
    def __init__(self, id, metaId,  subs, pros, subCoeff, proCoeff, lb,  ub):
        self.id = id
        self.metaId  = metaId
        self.substrates = set(subs)
        self.products = set(pros)
        self.subCoeff = subCoeff.copy()
        self.proCoeff = proCoeff.copy()
        self.lb = lb
        self.ub = ub
        self.overReaction = ''

    def getId(self):
        return self.id
        
    def getMetaId(self):
        return self.metaId
        
    def getlb(self):
        return self.lb;

    def getub(self):
        return self.ub;
        
    def getOverReaction(self):
        return self.overReaction;
    
    def getSubstrates(self):
        return self.substrates

    def getProducts(self):
        return self.products

    def getSubstrateCoeff(self):
        return self.subCoeff

    def getProductCoeff(self):
        return self.proCoeff

    def containsMetabolites(self, metabolites):
        for metabolite in metabolites:
            metId = metabolite.getId()
            for substrate in self.substrates:
                if substrate in metId:
                    return True
            for product in self.products:
                if product in metId:
                    return True
        return False

    def containsAllMetabolites(self, metabolites):
        for metabolite in metabolites:
            isThere = False
            allMetabolites = set(list(self.substrates) + list(self.products))
            for myMetabolite in allMetabolites:
                if myMetabolite in metabolite.getId():
                    isThere = True
            if not isThere:
                return False

        return True

    def isInModel(self, reactions):
        for reaction in reactions:
            if reaction.getId() == self.getId():
                self.overReaction = reaction.getId()
                return True
            sizeR1 = len(self.substrates) + len(self.products)
            sizeR2 = reaction.getNumReactants() + reaction.getNumProducts()
            if sizeR1 == sizeR2:
                if self.containsAllMetabolites(reaction.getListOfReactants())\
                and self.containsAllMetabolites(reaction.getListOfProducts()):
                    self.overReaction = reaction.getId()
                    return True

        return False
        
###############################################################################

class ModifyModel:
    def parseReactionBag(self, reactionBagPath):
            f = open(reactionBagPath)
            for s in f:
                if s.startswith("#"):
                    continue                
                key,name,eq,lb,ub= s.split(",")
                compounds = eq.split(" <=> ")
                lhs = ''
                rhs = ''
                if len(compounds) == 2:
                    lhs = compounds[0]
                    rhs = compounds[1]
                else:
                    lhs = compounds[0]
                    
                subs = set()
                pros = set()
                subCoeff = {}
                proCoeff = {}
                lhs, rhs = lhs.split(" + "), rhs.split(" + ")
                for mol in lhs:
                    vals = mol.split()
                    if len(vals) == 2:
                        coeff, mol = vals
                    elif len(vals) == 1:
                        coeff, mol = 1, vals[0]
                    else:
                        continue
                    subs.add(mol)
                    subCoeff[mol] = coeff
                for mol in rhs:
                    vals = mol.split()
                    if len(vals) == 2:
                        coeff, mol = vals
                    elif len(vals) == 1:
                        coeff, mol = 1, vals[0]
                    else:
                        continue
                    pros.add(mol)
                    proCoeff[mol] = coeff
                #print "creating reaction: "+ key    
                reaction = Reaction(key,  name, subs, pros, subCoeff,  proCoeff, lb, ub)
                #print "saving reaction: "+ key 
                self.reactions.append(reaction)
                #print "done" 
                    
    def parseBiomass(self,  biomassPath): 
        f = open(biomassPath)
        subs = set()
        pros = set()
        subCoeff = {}
        proCoeff = {}
        key = "biomasspseudoreaction"
        name = "Biomass_Pseudoreaction"
        lb = 0
        ub = 1000
        
        for s in f: 
            if s.startswith("#"):
                continue
            compound, coefficient= s.split(',')           
            if(float(coefficient) < 0):
                subs.add(compound)
                subCoeff[compound]=(abs(float(coefficient)))
            else:
                pros.add(compound)
                proCoeff[compound] = (abs(float(coefficient)))
        
        reaction = Reaction(key,  name, subs,  pros, subCoeff,  proCoeff, lb, ub)
        return reaction
        
        
    def parseBounds(self,  boundsPath,  model): 
        f = open(boundsPath) 
        for s in f: 
            if s.startswith("#"):
                continue
            id,  lb,  ub = s.split(',')
            reaction = model.getReaction(id)
            if reaction is not None:
                kinLaw = reaction.createKineticLaw()
                lbParameter= kinLaw.createLocalParameter()
                lbParameter.setId("LOWER_BOUND")
                lbParameter.setValue(float(lb))
                ubParameter= kinLaw.createLocalParameter()
                ubParameter.setId("UPPER_BOUND")
                ubParameter.setValue(float(ub))
                

    def addReactionInToModel(self, reaction, model):
            # Adds the metabolites that are not already in the model
            for metabolite in set(list(reaction.getSubstrates()) +
            list(reaction.getProducts())):
                if model.getSpecies(metabolite) is None:
                    sp = model.createSpecies()
                    sp.setName(metabolite)
                    sp.setId(metabolite)                   
                    sp.setCompartment('cytosol')
                    sp.setConstant(True)
            # Adds the reaction into the model
            r = model.createReaction()
            r.setId(reaction.getId())
            r.setMetaId("meta_"+ reaction.getMetaId())
            #r.setId(reaction.getId().replace("-", ""))

            # Adds reactants
            coeff = reaction.getSubstrateCoeff()
            for sp in reaction.getSubstrates():
                reactant = r.createReactant()
                reactant.setSpecies(sp)
                reactant.setStoichiometry(float(coeff[sp]))
                reactant.setConstant(True)
            # Adds products
            coeff = reaction.getProductCoeff()
            for sp in reaction.getProducts():
                product = r.createProduct()
                product.setSpecies(sp)
                product.setStoichiometry(float(coeff[sp]))
                product.setConstant(True)
                
            # Adds bounds
            lb = reaction.getlb()
            ub = reaction.getub()
            kinLaw = r.createKineticLaw()
            lbParameter= kinLaw.createLocalParameter()
            lbParameter.setId("LOWER_BOUND")
            lbParameter.setValue(float(lb))
            ubParameter= kinLaw.createLocalParameter()
            ubParameter.setId("UPPER_BOUND")
            ubParameter.setValue(float(ub))

            

   

    def __init__(self,  modelPath,  reactionPath,  biomassPath,  boundsPath,  outputPath):     
        self.reactions = []  
        self.model =''
        
        if(modelPath):
            #read model
            reader = libsbml.SBMLReader()
            document = reader.readSBML(modelPath)
            self.model = document.getModel()
            newDocument = document.clone();
            
             
            if(reactionPath):            
                self.parseReactionBag(reactionPath)
                if self.model == None:
                    return
                listOfReactions = self.model.getListOfReactions()
                for reaction in self.reactions:
                    if not reaction.isInModel(listOfReactions):
                        print "adding Reaction:" + reaction.getId()
                    else:
                        print "the reaction " + reaction.getId() + " is already on the model. OverWriting..."
                        self.model.getReaction(reaction.getOverReaction()).removeFromParentAndDelete()
                    self.addReactionInToModel(reaction, self.model)
                        
            if(biomassPath):
                print "adding biomass Reaction"
                reaction = self.parseBiomass(biomassPath)
                if not reaction.isInModel(self.model.getListOfReactions()):
                    self.addReactionInToModel(reaction, self.model)
                        
            if(boundsPath):           
                print "adding new bounds"
                self.parseBounds(boundsPath,  self.model)

            if(outputPath):
                print "Writing file"                
                newDocument.setModel(self.model)
                libsbml.writeSBMLToFile(newDocument, outputPath)   
                
            
            

if __name__ == "__main__":
    import sys,  getopt
    
    myopts, args = getopt.getopt(sys.argv[1:],"m:r:b:k:o:")
    modelPath = ''
    reactionPath = ''
    biomassPath = ''
    boundsPath = ''
    outputPath = ''
    model = ''  

    for o,  a in myopts:
        if o == '-m':
            modelPath = a
        elif o == '-r':
            reactionPath = a
        elif o =='-b':
            biomassPath = a
        elif o == '-k':
            boundsPath = a
        elif o == '-o':
            outputPath = a
        else:
            print("Usage: %s -m path_of_the_models -r path_of_the_reactions -br path_of_the_biomass_file -b path_of_the_bounds_file -o path_of_the_output_folder" % sys.argv[0])    
          
    ModifyModel(modelPath,  reactionPath,  biomassPath,  boundsPath,  outputPath)
