# kegg_ligand_parser.py - parse KEGG LIGAND reactions
#
# $Log: kegg_ligand_parser.py,v $
# Revision 1.5  2006/05/30 17:17:47  epitkane
# Added support for pubchem ids.
#
# Revision 1.4  2006/01/21 11:56:49  epitkane
# Added support for reaction enzymes, updated to reflect vocabulary.py
#
# Revision 1.3  2005/12/12 10:32:47  epitkane
# Added parseCompound()
#
# Revision 1.2  2005/10/19 10:47:12  epitkane
# ..
#
# Revision 1.1  2005/09/14 15:01:32  epitkane
# Initial version in cvs
#

import sys, re
from metabolism.network import *
from metabolism.vocabulary import *

class AtomGraph:
    def __init__(self, edges, carbonCount, pairToReaction, pairToSub, pairToPro, atomMap,
                 atomCount, atomTypes, atomMapTypes):
        self.edges = edges
        self.carbonCount = carbonCount
        self.pairToReaction = pairToReaction
        self.pairToSub = pairToSub              # pair id -> substrate of pair
        self.pairToPro = pairToPro
        self.atomMap = atomMap                  # pair id -> source atom -> target atom
        self.atomMapTypes = atomMapTypes        # pair id -> position -> atom type
        self.atomCount = atomCount              # compound id -> last atom index seen
        self.atomTypes = atomTypes              # compound id -> list of atom types

def parseAtomGraphFromRpair(f, forbiddenReactions = set(), verbose = False, errfile = None):
    """Parse an atom graph from rpair file (stream f)"""

    reAtomType = re.compile("(\D+)\d*\D*")

    g = {}
    pairToReaction = {}
    pairToSub = {}        # rpair -> substrate of rpair
    pairToPro = {}
    cost = {}
    atommap = {}
    numatoms = {}
    atomtypes = {}
    invalid = 0
    atommaptypes = {}

    numC = 0
    id = cmp1 = cmp2 = amap = None
    reactionBlock = 0
    res = None

    while 1:
        s = f.readline()
        #print s.strip("\n")

        if len(s) == 0:
            break
        
        if s[0] != " ":
            reactionBlock = 0
        
        if s.startswith("ENTRY "):
            id = s[12:20].strip()
            #print "ID===", id
        elif s.startswith("NAME"):
            cmp1, cmp2 = s[12:].strip().split("_")
            if cmp1 not in numatoms:
                numatoms[cmp1] = 0
                atomtypes[cmp1] = {}
            if cmp2 not in numatoms:
                numatoms[cmp2] = 0
                atomtypes[cmp2] = {}
        elif s.startswith("REACTION"):
            reactionBlock = 1
            res = set(s[12:].strip().split(" "))
            forbiddenResInRes = res.intersection(forbiddenReactions)
            if len(forbiddenResInRes) > 0:
                if verbose:
                    print "Warning: skipping rpair %s because of forbidden reaction(s) %s" % (id, ",".join(list(forbiddenResInRes)))
                if errfile:
                    errfile.write("%s\tForbidden reaction %s\n" % (id, ",".join(list(forbiddenResInRes))))
                invalid = 1
        elif s.startswith("ALIGN"):
            nalign = int(s.split()[1].strip())
            amap = {}
            amtypes = []
            for i in range(nalign):
                s = f.readline()

                #if cmp1 == "C00041" or cmp2 == "C00041":
                #    print id, cmp1, cmp2, s.strip()
                ###print s.strip()
                vals = s.strip().split()
                ix, first, second = vals[0:3]
                
                firstIx, firstType = first.split(":")
                
                m = reAtomType.search(firstType)
                #print firstType
                firstType = m.group(1)
                amtypes.append(firstType)

                firstIx = int(firstIx)
                if firstIx > numatoms[cmp1]:
                    numatoms[cmp1] = firstIx
                if firstType == "C":
                    numC += 1
                secondIx, secondType = second.split(":")
                m = reAtomType.search(secondType)
                if (m.group(1) != firstType):
                    if verbose:
                        print "Warning: invalid rpair: %s " % (id)
                    if errfile:
                        errfile.write("%s\tatom type mismatch\n" % (id))
                    invalid = 1
                    break
                    #print firstType, m.group(1), first, second, id
                    #sys.exit(1)
                secondIx = int(secondIx)

                # check against previous rpairs for these cmps
                if cmp1 in atomtypes and firstIx in atomtypes[cmp1]:
                    if atomtypes[cmp1][firstIx] != firstType:
                        if verbose:
                            print "Warning: mismatch in atom types: %s %d %s %s %s" % (cmp1, firstIx, id, atomtypes[cmp1][firstIx], firstType)
                        if errfile:
                            errfile.write("%s\tatom type mismatch\n" % (id))
                        invalid = 1
                if cmp2 in atomtypes and secondIx in atomtypes[cmp2]:
                    if atomtypes[cmp2][secondIx] != firstType:
                        if verbose:
                            print "Warning: mismatch in atom types: %s %d %s %s %s" % (cmp2, secondIx, id, atomtypes[cmp2][secondIx], secondType)
                        if errfile:
                            errfile.write("%s\tatom type mismatch\n" % (id))
                        invalid = 1
                        
                if invalid:
                    break

                if secondIx > numatoms[cmp2]:
                    numatoms[cmp2] = secondIx
                amap[int(firstIx)] = int(secondIx)

                atomtypes[cmp1][firstIx] = atomtypes[cmp2][secondIx] = firstType

        elif s[0] == " " and reactionBlock:
            res.update(set(s[12:].strip().split(" ")))
        elif s.startswith("///"):
            assert(id)
            assert(cmp1)
            assert(cmp2)

            if not invalid:
                if cmp1 not in g:
                    g[cmp1] = {}
                if cmp2 not in g[cmp1]:
                    g[cmp1][cmp2] = set()
                g[cmp1][cmp2].add(id)
                cost[id] = numC
                pairToReaction[id] = res
                pairToSub[id] = cmp1
                pairToPro[id] = cmp2
                atommap[id] = amap
                atommaptypes[id] = amtypes

            invalid = 0

            id = cmp1 = cmp2 = res = amap = None
            numC = 0
            
    return AtomGraph(g, cost, pairToReaction, pairToSub, pairToPro, atommap, numatoms, atomtypes, atommaptypes)

def checkReactionVsRpair(ag):
    reToPair = {}
    for pair in ag.pairToReaction:
        reactions = ag.pairToReaction[pair]
        if reactions:
            for re in reactions:
                if re not in reToPair:
                    reToPair[re] = set()
                reToPair[re].add(pair)

    pairToCmp = {}
    for cmp1 in ag.edges:
        for cmp2 in ag.edges[cmp1]:
            for pair in ag.edges[cmp1][cmp2]:
                if pair not in pairToCmp:
                    pairToCmp[pair] = {}
                pairToCmp[pair] = [cmp1, cmp2]

    for re in reToPair:
        good = True
        carbonsGood = True
        pairs = reToPair[re]
        mapped = {}
        for pair in pairs:
            #print pair
            #print pairToCmp[pair]
            src, tgt = pairToCmp[pair]
            for atom in ag.atomMap[pair]:
                satom = "%s-%d" % (src, atom)
                tatom = "%s-%d" % (tgt, ag.atomMap[pair][atom])
                if satom not in mapped:
                    mapped[satom] = 1
                else:
                    mapped[satom] += 1
                    good = False
                    if ag.atomTypes[src][atom] == "C":
                        carbonsGood = False
                if tatom not in mapped:
                    mapped[tatom] = 1
                else:
                    mapped[tatom] += 1
                    good = False
                    if ag.atomTypes[tgt][ag.atomMap[pair][atom]] == "C":
                        carbonsGood = False
        print "%s\t%s\t%s" % (re, good, carbonsGood)        
            
def parse_vocabulary_from_compound(f, callback = None):
    v = Vocabulary()
    count = 0
    while 1:
        count += 1

        if callback and (count % 1000) == 0:
            callback()
        
        block = parse_kegg_data_block(f)
        if block == None:
            break

        entry = block["ENTRY"].split()[0].strip()

        if "NAME" in block:
            nblock = block["NAME"].split("\n")
            first = 1
            for name in nblock:
                v.addSynonym(entry, name.strip(";").strip(), first)
                first = 0
        
    return v        

def parse_compound(f, callback = None):
    compounds = {}
    count = 0
    while 1:
        count += 1

        if callback and (count % 1000) == 0:
            callback()
        
        block = parse_kegg_data_block(f)
        if block == None:
            break

        entry = block["ENTRY"].split()[0].strip()
        if "NAME" in block:
            names = block["NAME"].split("\n")
            name = names[0].strip()
        else:
            name = ""

        if "FORMULA" in block:
            formula = block["FORMULA"].strip()
        else:
            formula = ""

        pubchemid = None
        if "DBLINKS" in block:
            dblinks = block["DBLINKS"].split("\n")
            for s in dblinks:
                m = re.search("PubChem: (\d+)", s)
                if m:
                    pubchemid = m.group(1)
                    #print entry, pubchemid

        consumers = []
        producers = []

        mol = Molecule(entry, consumers, producers)
        mol.names = names
        mol.formula = formula
        mol.pubchemid = pubchemid
        compounds[entry] = mol
        #print entry, compounds[entry].pubchemid
    return compounds

def parse_atoms_from_compound(f):

    num_comp = 0
    num_atominfo = 0

    atomtypes = {}
    names = {}

    sys.stderr.write("Atom types: ")
    sys.stderr.flush()
    
    count = 0
    while 1:
        count += 1

        if (count % 1000) == 0:
            sys.stderr.write(".")
            sys.stderr.flush()
        
        block = parse_kegg_data_block(f)
        if block == None:
            break

        num_comp += 1

        entry = block["ENTRY"].split()[0]

        if entry in atomtypes:
            print "Duplicate entry in atomtypes: %s" % (entry)
            return

        names[entry] = []
        if "NAME" in block:
            nblock = block["NAME"].split("\n")
            for name in nblock:
                names[entry].append(name.replace(";", "").strip())
            #print "%s\t%s" % (entry, names[entry])

        atomtypes[entry] = {}
        
        #print entry
        if "ATOM" in block:
            atoms = block["ATOM"]
            num_atominfo += 1

            atomlines = atoms.split("\n")
            for line in atomlines[1:]:
                elements = line.split()
                seq, type = elements[0], elements[2]
                atomtypes[entry][int(seq)] = type
           
    sys.stderr.write("\n")
    sys.stderr.write("%d compounds, %d with atom information\n" % (num_comp, num_atominfo))

    return atomtypes, names

# input: kegg ligand 'rpair', metabolic network
def parse_rpair(f, net):

    num_unknown_maps = 0
    count = 0
    
    sys.stderr.write("Atom maps: ")
    sys.stderr.flush()

    while 1:
        count += 1
        if (count % 1000) == 0:
            sys.stderr.write(".")
            sys.stderr.flush()
        
        rpblock = parse_kegg_data_block(f)
        if rpblock == None:
            break

        entry = rpblock["ENTRY"].split()[0]
        #print "Rpair:", entry
        compounds = rpblock["COMPOUND"].replace("\n", " ").split(" ")
        reactions = rpblock["REACTION"].replace("\n", " ").split(" ")
        ablock = rpblock["ALIGN"].split("\n")
        num_align = int(ablock[0])
        map = {}

        for line in ablock[1:]:
            elem = line.split()

            ix = int(elem[0])
            atom1, atom1type = elem[1].split(":")
            atom2, atom2type = elem[2].split(":")

            #print ix, atom1, atom1type, atom2, atom2type

            if atom1 == '' or atom2 == '' or atom1 == '0' or atom2 == '0': # incomplete/unknown mapping in KEGG, skip
                map = None
                break
            map["%s:%d" % (compounds[0], int(atom1))] = "%s:%d" % (compounds[1], int(atom2))

        # add this atom map to related reactions
        for rname in reactions:
            if rname == '':
                continue
            r = net.getReactionByName(rname)
            if r == None:
                print "Warning: can't find reaction %s for atom map addition - quit" % (rname)
                print rpblock
                sys.exit(1)
            else:
                if map != None:
                    r.addAtomMap(map)
                else:
                    r.atommap = None
                    r.unknownAtomMap = 1
                    num_unknown_maps += 1

        #print compounds
        #print reactions
        #print num_align

        #print "MAP:",map

    sys.stderr.write("\n")
    sys.stderr.write("%d atom maps in total, %d unknown/partial maps\n" % (count, num_unknown_maps))

def parse_network_dir(dir, progressCallback = None, addReverses = True):
    f = open("%s/reaction" % (dir))
    return parse_network(f, progressCallback, addReverses)

# input: kegg ligand 'reaction'
def parse_network(f, progressCallback = None, addReverses = True):
    net = MetabolicNetwork()
    r = parse_reaction(f)
    i = 0
    while r != None:
        net.addReaction(r)
        if addReverses:
            net.addReaction(r.reverse())
        r = parse_reaction(f)
        i += 1
        if progressCallback and (i % 1000) == 0:
            progressCallback()
    return net

# input: kegg ligand 'reaction'
def parse_reaction(f):
    data = parse_kegg_data_block(f)
    if data == None:
        return None
    
    if "ENTRY" not in data:
        raise Exception, "Expecting entry"

    entry = data["ENTRY"].split()[0]
    if "NAME" in data:
        name = data["NAME"].replace("\n", " ")
    else:
        name = "?"

    refid = -1
    bidir = 1
    pathname = "?"

    subs, prods = parse_reactants(data["EQUATION"].replace("\n", " "))

    if "ENZYME" in data:
        enzymes = data["ENZYME"].replace("\n", " ").split()
    else:
        enzymes = []

    r = Reaction(entry, name, refid, pathname, bidir, subs, prods, enzymes)

    # Note: atom mappings extracted from 'rpair' file
    #if "RPAIR" in data:
    #    r.rpair = data["RPAIR"].split(" ")
    #else:
    #    r.rpair = []

    #print "Reaction: %s: %s => %s" % (entry, subs, prods)

    return r

def parse_enzyme_to_gene(f):
    map = {}

    sys.stderr.write("Enzymes: ")
    sys.stderr.flush()
    
    e = parse_enzyme(f)

    map.update(e)
    
    i = 0
    while e != None:
        e = parse_enzyme(f)
        if e != None:
            map.update(e)
        i += 1
        if (i % 1000) == 0:
            sys.stderr.write(".")
            sys.stderr.flush()
    sys.stderr.write("\n")

    return map
    
def parse_enzyme(f):
    data = parse_kegg_data_block(f)
    if data == None:
        return None
    
    if "ENTRY" not in data:
        raise Exception, "Expecting entry"

    if "GENES" not in data:
        return {}

    id = data["ENTRY"][3:] # remove 'EC '

    map = {}
    map[id] = {}

    lines = data["GENES"].split("\n")
    organism = None
    for line in lines:
        if line.find(":") != -1:
            organism, genes = line.split(": ")
        else:
            assert(organism)
            genes = line

        genelist = genes.split()
        #print "%s : %s -> %s" % (id, organism, genelist)

        map[id][organism.lower()] = genelist

    return map

def parse_kegg_data_block(f):
    data = {}
    key = None
    while 1:
        stuff = parse_keyvalue(f)
        if stuff == None:         # key without value
            continue
        elif len(stuff) == 2:     # block starts
            key, value = stuff
            if key == None: # eof
                return None
            if key == "///":
                break
            data[key] = value.strip()
        elif len(stuff) == 1:   # block continues
            data[key] += "\n%s" % (stuff[0].strip())
    return data    

def parse_reactants(equation):
    subs, prods = [], []
    lhs, rhs = equation.split(" <=> ")
    subs = get_reactant_list(lhs)
    prods = get_reactant_list(rhs)
    return subs, prods
    
def get_reactant_list(stuff):
    list = []
    s = stuff.strip().split(" + ")
    #print "ReactantList:", s
    for sub in s:
        items = sub.strip().split(" ")
        if len(items) == 1:
            coeff = 1
            name = items[0].strip()
        else:
            try:
                coeff = int(items[0].strip())
            except ValueError, e:
                coeff = items[0].strip()
            name = items[1].strip()

        ix = name.find("(")
        if ix != -1:
            #print stuff
            #print name, coeff, ix
            coeff = name[ix + 1:-1]
            name = name[0:ix]
            #print name, coeff
           
        reactant = Reactant(name, coeff)
        list.append(reactant)
    return list
    
def parse_keyvalue(f):
    tries = 0
    while 1:
        s = f.readline()
        if s == "": # EOF
            return None, None
        s = s.strip("\n")
        if s != "":
            break
        tries += 1
        if tries > 100:
            print "Parse error: more than 100 empty lines.\nCheck input."
            sys.exit(2)
    #if len(s) == 0:
    #    return None, None
    if s.startswith("///"):
        return "///", None
    # starts block?
    #               key     value 
    m = re.search("^(\w+)\W+(.+)$", s)
    if m == None:
        # continues block or is key without value
        m = re.search("^\s+(\S.*)$", s)
        if m == None:
            # key without value
            return None
        # block continues
        return [m.group(1)]
    return m.group(1), m.group(2)

def parseReactionList(f):
    block = 0
    rl = {}
    for s in f:
        if block:
            if s[0] == " ":
                eqn += s.rstrip()[11:]
                continue
            else:
                block = 0
        if s.startswith("ENTRY"):
            id = s.split()[1]
        elif s.startswith("EQUATION"):
            eqn = s.strip()[12:]
            block = 1
        elif s.startswith("///"):
            rl[id] = eqn
    return rl

def kegg_ligand_parser_test():
    f = open(sys.argv[1])
    compounds = parse_compound(f, None)
    f.close()

def printReactionList(fn):
    f = open(fn)
    rl = parseReactionList(f)
    for re in rl:
        sys.stdout.write("%s\t%s\n" % (re, rl[re]))

if __name__ == '__main__':
    #ag = parseAtomGraphFromRpair(open("/group/home/icomic/data/kegg/ligand/LATEST/rpair"))
    #checkReactionVsRpair(ag)
    #kegg_ligand_parser_test()
    f = open("/group/home/icomic/data/kegg/ligand/LATEST/reaction")
    rl = parseReactionList(f)
    for re in rl:
        mult = "Single"
        lhs, rhs = rl[re].split(" <=> ")
        print lhs, rhs
        reacts = lhs.split(" + ")
        reacts.extend(rhs.split(" + "))
        for r in reacts:
            vals = r.split(" ")
            if len(vals) > 1:
                mult = "Multiple"
        sys.stdout.write("%s\t%s\t%s\n" % (re, rl[re], mult))
