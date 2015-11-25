#
# network.py - stoichiometric metabolic network and related classes
#
# Package: metabolism
# Provides:
#    Molecule
#    Reactant
#    Reaction
#    MetabolicNetwork
#
# $Log: network.py,v $
# Revision 1.4  2006/01/21 11:53:50  epitkane
# Added support for reaction enzymes
#
# Revision 1.3  2005/12/12 10:33:47  epitkane
# Some minor changes
#
# Revision 1.2  2005/10/19 10:47:04  epitkane
# ..
#
# Revision 1.1  2005/09/14 15:01:21  epitkane
# Initial version in cvs
#

import sys
#from metabolism.utils.altsets import Set
#from sets import Set

class Molecule:

    def __init__(self, name, consumers = [], producers = []):
        assert(name)
        self.name = name

        self.consumers = {}
        for c in consumers:
            self.consumers[c.name] = c

        self.producers = {}
        for p in producers:
            self.producers[p.name] = c

    def __cmp__(self, other):
        if not isinstance(other, Molecule):
            return 1
        return cmp(self.name, other.name)

    def __hash__(self):
        assert(self.name)
        return hash(self.name)

    def __str__(self):
        assert(self.name)
        return self.name

    def __repr__(self):
        assert(self.name)
        return self.name

class Reactant:

    def __init__(self, name, coeff = 1):
        assert(name)
        self.name = name
        self.coeff = coeff

    def copy(self):
        r = Reactant(self.name, self.coeff)
        return r

    def __cmp__(self, other):
        return cmp(self.name, other.name)

    def __str__(self):
        return "%s*%s" % (self.coeff, self.name)

    def __repr__(self):
        return "%s*%s" % (self.coeff, self.name)

    def __hash__(self):
        return hash(self.name)

class Reaction:

    def __init__(self, name, description, refid, pathname, bidirectional, subs, prods, enzymes = []):
        self.name = name
        self.description = description
        self.refid = refid
        self.pathname = pathname
        self.bidir = bidirectional

        self.substrates = {}
        for s in subs:
            self.substrates[s.name] = s

        self.products = {}
        for p in prods:
            self.products[p.name] = p

        self.atommap = {}
        self.unknownAtomMap = 0 # 0=still possible to construct full map, 1=full map impossible to construct

        self.enzymes = list(enzymes)

    def reverse(self):
        r = Reaction("%s_rev" % (self.name), self.description, self.refid, self.pathname, self.bidir, [], [], [])

        for s in self.substrates:
            r.products[s] = self.substrates[s].copy()

        for p in self.products:
            r.substrates[p] = self.products[p].copy()

        r.atommap = self.atommap.copy()
        r.unknownAtomMap = self.unknownAtomMap

        r.enzymes = list(self.enzymes)

        return r
        

    def copy(self):
        r = Reaction(self.name, self.description, self.refid, self.pathname, self.bidir, [], [], [])

        for s in self.substrates:
            r.substrates[s] = self.substrates[s].copy()

        for p in self.products:
            r.products[p] = self.products[p].copy()

        r.atommap = self.atommap.copy()
        r.unknownAtomMap = self.unknownAtomMap

        r.enzymes = list(self.enzymes)

        return r

    def addAtomMap(self, map):
        "Add the given atom map to any existing map"

        if self.unknownAtomMap: # cannot construct complete map so skip this entry too
            return
        
        # map: compounds x atoms -> compounds x atoms
        # e.g., C00010:1 -> C01234:22
        # note: mapped metabolite can be either substrate or product,
        # reverse if product -> substrate
        for k in map:

            cmp, atom = k.split(":")
            if cmp not in self.substrates:
                # reverse to get substrate -> product
                key = map[k]
                value = k
                #print "Reversing map %s -> %s" % (k, map[k])
            else:
                # correct order
                key = k
                value = map[k]
                #print "Nice map %s -> %s" % (k, map[k])
            
            if key not in self.atommap:
                self.atommap[key] = []

            self.atommap[key].append(value)

    def getEquation(self):
        s = ""
        for m in self.substrates:
            s += "%s + " % (m)
        s = s.strip("+ ")
        s = "%s " % (s)
        if self.bidir:
            s += "<"
        s += "=>"
        for m in self.products:
            s += " %s +" % (m)
        s = s.strip(" +")
        return s

    def strWithReplacement(self, vocabulary):
        s = str(self.name) + ": "
        for m in self.substrates:
            name = vocabulary.getNonIdSynonym(m) or "<noname>"
            s += "%s (%s) " % (name, m)

        if self.bidir:
            s += "<"

        s += "=>"

        for m in self.products:
            name = vocabulary.getNonIdSynonym(m) or "<noname>"
            s += " %s (%s)" % (name, m)
        
        return s

    def __cmp__(self, other):
        if not isinstance(other, Reaction):
            return -1
        else:
            return cmp(self.name, other.name)

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        s = str(self.name) + ": "
        for m in self.substrates:
            s += "%s " % (m)
        if self.bidir:
            s += "<"
        s += "=>"
        for m in self.products:
            s += " %s" % (m)
        return s

    def __repr__(self):
        return self.name

class MetabolicNetwork:

    def __init__(self):
        self.id = None
        self.name = ""
        self.description = ""
        self.author = ""
        self.create_date = "" 
        self.organism = ""
        self.ispublic = 0
        self.reactions = {}
        self.molecules = {}
        self.catalysedBy = {}
        self.catalyses = {}
        self.enzymeORFs = {}
        self.geneORFs = {}

    def copy(self):
        net = MetabolicNetwork()
        net.id = self.id
        net.name = self.name
        net.description = self.description
        net.author = self.author
        net.create_date = self.create_date
        net.organism = self.organism
        net.ispublic = self.ispublic
        for r in self.reactions:
            net.addReaction(self.reactions[r].copy())
        net.catalysedBy = self.catalysedBy.copy()
        net.catalyses = self.catalyses.copy()
        net.enzymeORFs = self.enzymeORFs.copy()
        net.geneORFs = self.geneORFs.copy()

        # note: removed self.atommap.copy
        
        return net

    def merge(self, net):
        """Merge the given network into self. Do not add reactions already in self."""

        for r in net.reactions:
            if not self.hasReaction(r):
                self.addReaction(net.reactions[r])

    def addReverseReactions(self):
        adds = []
        for r in self.reactions:
            rev = self.reactions[r].reverse()
            if rev.name not in self.reactions:
                adds.append(rev)
        for a in adds:
            self.addReaction(a)

    def getSubset(self, R):
        net = MetabolicNetwork()
        for r in R:
            if r in self.reactions:
                net.addReaction(self.reactions[r].copy())
        return net

    def addCatalysis(self, enzyme, reaction):
        assert(enzyme) # enzyme not given
        assert(reaction) # reaction not given
        if enzyme not in self.catalyses:
            self.catalyses[enzyme] = [reaction]
        else:
            if reaction not in self.catalyses[enzyme]:
                self.catalyses[enzyme].append(reaction)
        if reaction not in self.catalysedBy:
            self.catalysedBy[reaction.name] = [enzyme]
        else:
            if enzyme not in self.catalysedBy[reaction.name]:
                self.catalysedBy[reaction.name].append(enzyme)            

    def addORFForEnzyme(self, orf, ec):
        if ec not in self.enzymeORFs:
            self.enzymeORFs[ec] = [orf]
        else:
            if orf not in self.enzymeORFs[ec]:
                self.enzymeORFs[ec].append(orf)            

    def addORFForGene(self, orf, gene):
        if gene not in self.geneORFs:
            self.geneORFs[gene] = [orf]
        else:
            if orf not in self.geneORFs[gene]:
                self.geneORFs[gene].append(orf)            

    def getReactionByName(self, rname):
        if rname in self.reactions:
            return self.reactions[rname]
        else:
            return None
        
    def getMoleculeByName(self, mname):
        if mname in self.molecules:
            return self.molecules[mname]
        else:
            return None

    def getReactionNames(self):
        A = []
        for r in self.reactions:
            A.append(r)
        return A

    def getMoleculeNames(self):
        A = []
        for m in self.molecules:
            A.append(m)
        return A

    def removeReaction(self, rname):
        reaction = self.getReactionByName(rname)
        if reaction == None:
            return

        # remove reactants; remove molecule iff not used by any reaction anymore
        for m in reaction.substrates:
            molecule = self.getMoleculeByName(m)
            assert(molecule)
            del molecule.consumers[rname]
            if len(molecule.consumers) == 0 and len(molecule.producers) == 0:
                del self.molecules[m]

        for m in reaction.products:
            molecule = self.getMoleculeByName(m)
            assert(molecule)
            del molecule.producers[rname]
            if len(molecule.consumers) == 0 and len(molecule.producers) == 0:
                del self.molecules[m]

        del self.reactions[rname]

    def removeReactions(self, R):
        for r in R:
            self.removeReaction(r)

    def removeMolecule(self, mname):
        molecule = self.getMoleculeByName(mname)

        if molecule == None:
            return

        # remove molecule from reaction reactant lists
        for rn in molecule.consumers:
            r = molecule.consumers[rn]
            del r.substrates[mname]

        for rn in molecule.producers:
            r = molecule.producers[rn]
            del r.products[mname]

        del self.molecules[mname]

    def removeMolecules(self, M):
        for m in M:
            self.removeMolecule(m)
            
    def hasMolecule(self, mname):
        return mname in self.molecules

    def hasReaction(self, rname):
        return rname in self.reactions

    def addReactions(self, R):
        for r in R:
            self.addReaction(r)

    def addReaction(self, r):
        assert(not self.hasReaction(r.name))

        self.reactions[r.name] = r

        for mname in r.substrates:
            m = r.substrates[mname]
            molecule = self.getMoleculeByName(m.name)
            if molecule == None:
                molecule = Molecule(m.name)
                self.molecules[m.name] = molecule
            molecule.consumers[r.name] = r

        for mname in r.products:
            m = r.products[mname]
            molecule = self.getMoleculeByName(m.name)
            if molecule == None:
                molecule = Molecule(m.name)
                self.molecules[m.name] = molecule
            molecule.producers[r.name] = r

    def getExternalMolecules(self):
        E = []
        for m in self.molecules:
            mol = self.molecules[m]
            if len(mol.consumers) == 0:
                type = 0 # product
            elif len(mol.producers) == 0:
                type = 1 # substrate
            else:
                continue
                #type = 2 # internal
            E.append((m, type))

        return E

    def size(self):
        return len(self.reactions)

    def printMolecules(self):
        for m in self.molecules:
            mol = self.molecules[m]
            sys.stdout.write(mol.name + "\t")
            
            for consumer in mol.consumers:
                sys.stdout.write(" " + consumer.name)
            sys.stdout.write("\t")
            for producer in mol.producers:
                sys.stdout.write(" " + producer.name)
            sys.stdout.write("\n")

    def catalysisToStr(self):
        s = ""
        for ec in self.catalyses:
            s += ec + " -> " + str(self.catalyses[ec]) + "\n"
        return s

    def ecOrfToStr(self):
        s = ""
        for ec in self.enzymeORFs:
            s += ec + " -> " + str(self.enzymeORFs[ec]) + "\n"
        return s

    def geneOrfToStr(self):
        s = ""
        for gene in self.geneORFs:
            s += gene + " -> " + str(self.geneORFs[gene]) + "\n"
        return s

    def __str__(self):
        s = ""
        for rname in self.reactions:
            r = self.reactions[rname]
            s += "%s\n" % (str(r))
        return s.strip()

    def __repr__(self):
        s = "MetabolicNetwork: " + str(len(self.reactions)) + " reactions, " + str(len(self.molecules)) + " molecules"
        return s
