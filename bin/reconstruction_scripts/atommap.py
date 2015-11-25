#!/usr/bin/env python

import sys

class Reaction:
    """Metabolic reaction represented as an atom mapping.

An atom mapping is a many-to-many mapping of source atoms to target atoms. Reactions with singular reactants have one-to-one mappings. However, a mapping involving multiple reactants of the same metabolite maps the same atom more than once."""

    def __init__(self, id, subs, pros, sub_coeff, pro_coeff, maps, rev_maps = None, spontaneous = False):
        self.id = id
        self.substrates = set(subs)
        self.products = set(pros)
        self.sub_coeff = sub_coeff.copy()
        self.pro_coeff = pro_coeff.copy()
        self.spontaneous = spontaneous
        self.maps = maps  # source atom -> set(target atoms)
        # rev_maps: target atom -> set(source atoms)
        if rev_maps == None:
            self.rev_maps = self.reverse_maps(maps)
        else:
            self.rev_maps = rev_maps

    def basename(self):
        return self.id.split("_", 1)[0]

    def reverse(self):
        """Return a copy of the reaction reversed."""
        rid = "%s_rev" % (self.id)
        return Reaction(rid, self.products, self.substrates, self.pro_coeff, self.sub_coeff, self.reverse_maps(self.maps), self.maps)

    def reverse_maps(self, maps):
        rmaps = {}
        for u in maps:
            for v in maps[u]:
                if v not in rmaps:
                    rmaps[v] = set()
                rmaps[v].add(u)
        return rmaps

    def map_hash(self):
        h = 0
        for u in self.maps:
            h2 = 0
            for v in self.maps[u]:
                h2 = hash((h2, v))
            h = hash((h, u, h2))
        return h

    # def remove_duplicate_mappings(self):
    #     """Remove mappings that are identical. Such mappings occur due to atom types being removed from data during parsing (e.g. oxygens and hydrogens). Retained mappings are renumbered from one."""
    #     hashes = set()
    #     new_maps = {}
    #     new_rev_maps = {}
    #     c = 1
    #     h = 0
    #     for u in self.maps:
    #         h2 = 0
    #         for v in self.maps[u]:
    #             h2 = hash((h2, v))
    #         h = hash((h, h2))
    #     if h in hashes:
    #         #print "Removing mapping %s" % (mid)
    #         pass
    #     else:
    #         hashes.add(h)
    #         new_maps[c] = self.maps
    #         new_rev_maps[c] = self.rev_maps
    #         c += 1
    #     # print "%s\t%d\t%d" % (self.id, len(self.maps), len(new_maps))
    #     # n = len(self.maps)
    #     # global num_maps, num_pruned_maps
    #     # if n not in num_maps:
    #     #     num_maps[n] = 1
    #     # else:
    #     #     num_maps[n] += 1
    #     # n = len(new_maps)
    #     # if n not in num_pruned_maps:
    #     #     num_pruned_maps[n] = 1
    #     # else:
    #     #     num_pruned_maps[n] += 1
    #     #if len(new_maps) > 1:
    #     #    for m in new_maps:
    #     #        print m, new_maps[m]
    #     #    import sys
    #     #    sys.exit()
    #     self.maps = new_maps
    #     self.rev_maps = new_rev_maps

    def __eq__(self, other):
        return self.id == self.other

    def __hash__(self, other):
        return hash(self.id)

    def description(self):
        return "%s: %s <=> %s" % (self.id, " + ".join(map(lambda x: " ".join(map(str, x)), self.sub_coeff.items())), " + ".join(map(lambda x: " ".join(map(str, x)), self.pro_coeff.items())))


    def __str__(self):
        s = "Reaction %s. Map: %s\n" % (self.id, self.maps)
        return s

    def __repr__(self):
        return "%s(%d,%d,%d)" % (self.id, len(self.substrates), len(self.products), len(self.maps))

class AtomMappingParser:
    """A parser for the CSBB atom mapping format."""

    def __init__(self):
        self.reactions = {}
        self.base_reactions = {}
        self.atom_index = {}        # "mol-atom" -> atom_index
        self.rev_atom_index = {}    # atom_index -> (mol, atom)
        self.atom_types = {}        # mol -> atom -> atom_type
        self.disregarded_atom_types = None
        self.accepted_atom_types = None
        self.next_atom_index = 1
        self.spontaneous_reactions = set()   # base_reaction -> boolean

    def delete_reactions(self, R):
        for r in R:
            if r in self.base_reactions:
                for r2 in self.base_reactions[r]:
                    del self.reactions[r2]
                del self.base_reactions[r]
                #print "Deleted %s" % (r)
            else:
                print "Cannot delete %s - not in stoichiometry" % (r)

    def get_atom_index(self, mol, atom):
        """(mol, atom) -> atom_index"""
        key = "%s-%s" % (mol, atom)
        if key not in self.atom_index:
            self.atom_index[key] = self.next_atom_index
            self.rev_atom_index[self.next_atom_index] = (mol, atom)
            self.next_atom_index += 1
        #print "map", mol, atom, self.atom_index[key]
        return self.atom_index[key]

    def get_mol_and_atom(self, atom_index):
        """atom_index -> (mol, atom)"""
        return self.rev_atom_index[atom_index]

    def has_mol_and_atom(self, atom_index):
        return atom_index in self.rev_atom_index

    def register_atom_type(self, mol, index, atype):
        if mol not in self.atom_types:
            self.atom_types[mol] = {}
        if index not in self.atom_types[mol]:
            self.atom_types[mol][index] = atype
        elif self.atom_types[mol][index] != atype:
            raise Exception("Atom type mismatch: %s-%s %s vs %s" % (mol, index, self.atom_types[mol][index], atype))       

    def parse_mapping(self, f, index2mol, duplicate_reactants):
        s = f.readline().strip()
        if s == "":
            return None, None
        key, map_id = s.split()
        if key != "MAPPING":
            raise Exception("Parse error: expecting MAPPING, got %s" % (key))
        map_id = int(map_id)
        f.readline() # ATOM ROW MOL ATOM ROW MOL
        mappings = {}
        while True:
            s = f.readline()
            if s == None:
                break
            if s.startswith("BONDS"):
                return map_id, mappings
            #print s.strip()
            atom1, row1, mol1, atom2, row2, mol2 = s.strip().split()
            if atom1 != atom2:
                raise Exception("Atom type mismatch: %s" % (s))
            row1, mol1, row2, mol2 = int(row1), int(mol1), int(row2), int(mol2)
            mol1, mol2 = index2mol[mol1], index2mol[mol2]
            if mol1 not in self.atom_types:
                self.atom_types[mol1] = {}
            if mol2 not in self.atom_types:
                self.atom_types[mol2] = {}
            if self.disregarded_atom_types != None and atom1 in self.disregarded_atom_types:
                #print "Skipping atom type: ", atom1, mol1, row1
                continue
            if self.accepted_atom_types != None and atom1 not in self.accepted_atom_types:
                continue
            self.register_atom_type(mol1, row1, atom1)
            self.register_atom_type(mol2, row2, atom2)
            atom_ix1 = self.get_atom_index(mol1, row1)
            atom_ix2 = self.get_atom_index(mol2, row2)
            if duplicate_reactants == False and atom_ix1 in mappings:
                raise Exception("Duplicate mapping %s -> %s" % (self.get_mol_and_atom(atom_ix1), self.get_mol_and_atom(atom_ix2)))
            if atom_ix1 not in mappings:
                mappings[atom_ix1] = set()
            mappings[atom_ix1].add(atom_ix2)
        raise Exception("Parse error: expecting BONDS")

    def parse_reaction_graph(self, f):
        f.readline() # ATOM1 MOL ATOM2 MOL ORDER TYPE
        graph = {}
        while True:
            s = f.readline().strip()
            if s == "":
                return graph
            atom1, mol1, atom2, mol2, order, btype = s.split()
        raise Exception("Parse error: expecting empty line after reaction graph")

    def parse_mapping_and_reaction_graph(self, f, index2mol, duplicate_reactants):
        maps = {}
        while 1:
            map_id, mapping = self.parse_mapping(f, index2mol, duplicate_reactants)
            if map_id == None:
                #print "MAPS", len(maps)
                return maps
            reaction_graph = self.parse_reaction_graph(f)
            #if duplicate_reactants:
            #    maps[map_id] = {}
            #else:
            maps[map_id] = mapping

    def parse_reaction(self, f):
        key, rid = f.readline().strip().split()
        if key != "REACTION":
            raise Exception("Parse error: expecting REACTION, got %s" % (key))
        key, eq = f.readline().strip().split(None, 1)
        if key != "EQUATION":
            raise Exception("Parse error: expecting EQUATION, got %s" % (key))
        key, indices = f.readline().strip().split(None, 1)
        if key != "INDICES":
            raise Exception("Parse error: expecting 'INDICES', got '%s'" % (key))
        key, cost = f.readline().strip().split(None, 1)
        if key != "COST":
            raise Exception("Parse error: expecting 'COST', got '%s'" % (key))
        key, time = f.readline().strip().split(None, 1)
        if key != "TIME":
            raise Exception("Parse error: expecting 'TIME', got '%s'" % (key))

        lhs, rhs = indices.split("=>")
        lhs, rhs = lhs.split("+"), rhs.split("+")

        leq, req = eq.split("=>")
        leq, req = leq.split("+"), req.split("+")
        index2mol = {}
        subs = set()
        pros = set()
        sub_coeff = {}
        pro_coeff = {}
        duplicate_reactants = False
        for i, mol in enumerate(leq):
            mol = mol.replace("\t", "").replace(" ", "")
            index = int(lhs[i].replace("\t", "").replace(" ", ""))
            if mol in subs:
                duplicate_reactants = True
            subs.add(mol)
            if mol not in sub_coeff:
                sub_coeff[mol] = 1
            else:
                sub_coeff[mol] += 1
            index2mol[index] = mol
        for i, mol in enumerate(req):
            mol = mol.replace("\t", "").replace(" ", "")
            index = int(rhs[i].replace("\t", "").replace(" ", ""))
            if mol in pros:
                duplicate_reactants = True
            pros.add(mol)
            if mol not in pro_coeff:
                pro_coeff[mol] = 1
            else:
                pro_coeff[mol] += 1
            index2mol[index] = mol

        f.readline()

        maps = self.parse_mapping_and_reaction_graph(f, index2mol, duplicate_reactants)
        reactions = []
        mapkeys = maps.keys()
        mapkeys.sort()
        #for m in mapkeys[0:1]:
        for m in mapkeys:
            name = "%s_%s" % (rid, m)
            re = Reaction(name, subs, pros, sub_coeff, pro_coeff, maps[m])
            reactions.append(re)
        return reactions

    def parse_reaction_dir(self, ddir, 
                           accepted_atom_types = set(["C"]), 
                           skip_atom_types = None):
        import os
        fns = os.listdir(ddir)
        fns.sort()
        self.reactions = {}
        self.base_reactions = {}
        #self.disregarded_atom_types = set(["H"])
        #self.disregarded_atom_types = set(["O", "H"])
        #self.disregarded_atom_types = set(["O", "H", "P"])
        #self.disregarded_atom_types = set(["O", "H", "C", "N"])
        #self.disregarded_atom_types = set(["O", "H", "C", "R", "S", "Mn"])
        #self.disregarded_atom_types = set(["O", "H", "N", "P", "R", "S", "Mn"])
        self.disregarded_atom_types = skip_atom_types
        self.accepted_atom_types = accepted_atom_types
        global num_maps, num_pruned_maps
        num_maps = {}
        num_pruned_maps = {}
        for i, fn in enumerate(fns):
            f = open("%s/%s" % (ddir, fn))
            try:
                res = self.parse_reaction(f)
            except:
                print "Warning: cannot parse reaction %s" % (fn)
                continue

            map_hashes = set()
            for r in res:
                h = r.map_hash()
                if h in map_hashes:
                    #print "Removing duplicate mapping", r.id
                    continue
                map_hashes.add(h)
                rev = r.reverse()
                self.reactions[r.id] = r
                self.reactions[rev.id] = rev
                #print rev
                #print "->", r.id
                #print "->", rev.id

                base = r.basename()
                if base not in self.base_reactions:
                    self.base_reactions[base] = []
                self.base_reactions[base].append(r.id)
                self.base_reactions[base].append(rev.id)

                #print r
                #print self.atom_types

            #if i == 1000000:
            #    print "atommap.py: BREAK"
            #    break

        # k = num_maps.keys()
        # k.sort(lambda x, y: cmp(int(x), int(y)))
        # import sys
        # sys.stderr.write("#OrigMapLen Count\n")
        # for u in k:
        #     sys.stderr.write("%s\t%s\n" % (u, num_maps[u]))
        # k = num_pruned_maps.keys()
        # k.sort(lambda x, y: cmp(int(x), int(y)))
        # sys.stderr.write("#PrunedMapLen Count\n")
        # for u in k:
        #     sys.stderr.write("%s\t%s\n" % (u, num_pruned_maps[u]))

    def parse_ligand_reaction(self, ddir, augment_only = True):
        """Add reactions by parsing KEGG LIGAND reaction file.
        If augment_only == True, a reaction is added only if it is not
        found in base_reactions. Otherwise the reaction is always added.
        """
        f = open("%s/reaction" % (ddir))
        block_eq = 1
        spontaneous = False
    
        for s in f:
            if s.startswith("///"):
                #print "Entry (%s), eq (%s)" % (entry, eq)
                if augment_only and entry in self.base_reactions:
                    #print "Not adding %s" % (entry)
                    continue
                #else:
                #    print "%s" % (entry)

                lhs, rhs = eq.split(" <=> ")
                subs = set()
                pros = set()
                sub_coeff = {}
                pro_coeff = {}
                lhs, rhs = lhs.split(" + "), rhs.split(" + ")
                for mol in lhs:
                    vals = mol.split()
                    if len(vals) == 2:
                        coeff, mol = vals
                    else:
                        coeff, mol = 1, vals[0]
                    subs.add(mol)
                    sub_coeff[mol] = coeff
                    if mol not in self.atom_types:
                        self.atom_types[mol] = {}

                for mol in rhs:
                    vals = mol.split()
                    if len(vals) == 2:
                        coeff, mol = vals
                    else:
                        coeff, mol = 1, vals[0]
                    pros.add(mol)
                    pro_coeff[mol] = coeff
                    if mol not in self.atom_types:
                        self.atom_types[mol] = {}

                maps = {}

                name = "%s_0" % (entry)

                re = Reaction(name, subs, pros, sub_coeff, pro_coeff, maps, spontaneous = spontaneous)
                rev = re.reverse()

                self.reactions[re.id] = re
                self.reactions[rev.id] = rev
                base = re.basename()
                if base not in self.base_reactions:
                    self.base_reactions[base] = []
                self.base_reactions[base].append(re.id)
                self.base_reactions[base].append(rev.id)

                block_eq = 0
                spontaneous = False
                continue
            key, value = s[0:12].strip(), s[12:].strip()
            if key == "ENTRY":
                entry = value.split()[0]
            elif key == "EQUATION":
                block_eq = 1
                eq = value
            elif key == "COMMENT":
                if "spontaneous" in value or "non enzymatic" in value or "Non-enzymatic" in value or "non-enzymatic" in value:
                    self.spontaneous_reactions.add(entry)
                    spontaneous = True
            elif key != "":
                block_eq = 0
            elif block_eq:
                eq += " " + value

if __name__ == "__main__":
    import sys
    #f = open(sys.argv[1])
    #parse_reaction(f)
    amp = AtomMappingParser()
    amp.parse_reaction_dir(sys.argv[1])
    amp.parse_ligand_reaction(sys.argv[2])


    
