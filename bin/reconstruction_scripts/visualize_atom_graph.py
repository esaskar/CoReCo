#!/usr/bin/env python

import sys, tempfile, subprocess

sys.path.append("../model_training_scripts/")
import common

import atommap

DRAW_SEPARATELY_LIST = "mols-to-draw-separately.txt"

class AtomGraphVisualizer:

    def __init__(self, amp, cdir):
        """amp AtomMappingParser"""
        self.amp = amp
        self.separate_mols = common.read_set(open("%s/%s" % (cdir, DRAW_SEPARATELY_LIST)))
        self.sep_id = 1

    def call_dot(self, imagetype, ifn, ofn):
        subprocess.call("dot -T%s %s -o %s" % (imagetype, ifn, ofn), shell = True)

    def get_mol_draw_id(self, mid):
        s = "%s_%s" % (mid, self.sep_id)
        self.sep_id += 1
        return s

    def get_mol_def(self, mid, mol, mol_names, mcolors, mfillcolors):
        if mol in mol_names:
            mlabel = mol_names[mol]
        else:
            mlabel = mol
        if mol in mcolors:
            mcolor = mcolors[mol]
        else:
            mcolor = "black"
        if mol in mfillcolors:
            fcolor = mfillcolors[mol]
        else:
            fcolor = "white"
        return "   %s [label=\"%s\",color=\"%s\",fillcolor=\"%s\",style=\"filled\"];\n" % (mid, mlabel, mcolor, fcolor)

    def draw_metabolic_network(self, reactions, ofn, 
                               rcolors = {}, rfillcolors = {}, rlabels = {},
                               mcolors = {}, mfillcolors = {}, 
                               mol_names = {},
                               imagetype = "png", glabel = None):
        of = tempfile.NamedTemporaryFile()
        of.write("digraph atomgraph {\n")
        if glabel != None:
            of.write("  label=\"%s\";\n" % (glabel))

        mols = set()
        for r in reactions:
            re = self.amp.reactions[r]
            mols.update(re.substrates)
            mols.update(re.products)

        for mol in mols:
            if mol in self.separate_mols:
                continue
            of.write(self.get_mol_def(mol, mol, mol_names, mcolors, mfillcolors))
        for r in reactions:
            re = self.amp.reactions[r]
            if r in rcolors:
                color = rcolors[r]
            else:
                color = "black"
            if r in rfillcolors:
                fcolor = rfillcolors[r]
            else:
                fcolor = "black"
            if r in rlabels:
                label = rlabels[r]
            else:
                label = "%s" % (r)
            of.write("   %s [label=\"%s\",shape=box,style=\"filled\",color=\"%s\",fillcolor=\"%s\"];\n" % (r, label, color, fcolor))
            for sub in re.substrates:
                if sub in self.separate_mols:
                    mid = self.get_mol_draw_id(sub)
                    of.write(self.get_mol_def(mid, sub, mol_names, mcolors, mfillcolors))
                else:
                    mid = sub
                of.write("   %s->%s;\n" % (mid, r))
            for pro in re.products:
                if pro in self.separate_mols:
                    mid = self.get_mol_draw_id(pro)
                    of.write(self.get_mol_def(mid, pro, mol_names, mcolors, mfillcolors))
                else:
                    mid = pro
                of.write("   %s->%s;\n" % (r, mid))
        of.write("}\n")
        of.flush()
        self.call_dot(imagetype, of.name, "%s.%s" % (ofn, imagetype))
        of.close()

    def draw_atom_graph(self, reactions, atoms, ofn, 
                        rcolors = {}, acolors = {}, 
                        imagetype = "png", 
                        draw_edge_labels = False,
                        glabel = None,
                        atom_edge_colors = {}):
        of = tempfile.NamedTemporaryFile()
        of.write("digraph atomgraph {\n")
        if glabel != None:
            of.write("  label=\"%s\";\n" % (glabel))

        cid = 0
        #print "atoms", atoms
        nodes = set()
        molatoms = {}
        edges = []
        for r in reactions:
            re = self.amp.reactions[r]
            for u in re.maps:
                if u not in atoms:
                    continue
                mol1, atom1 = self.amp.get_mol_and_atom(u)
                if mol1 not in molatoms:
                    molatoms[mol1] = set()
                molatoms[mol1].add(u)
                nodes.add(u)
                for v in re.maps[u]:
                    if v not in atoms:
                        continue
                    nodes.add(v)
                    edges.append((u, v))
                    mol2, atom2 = self.amp.get_mol_and_atom(v)
                    if mol2 not in molatoms:
                        molatoms[mol2] = set()
                    molatoms[mol2].add(v)                   

        #atomshape = "circle"
        atomshape = "point"

        for mol in molatoms:
            of.write(" subgraph cluster_%d {\n   label=\"%s\";\n" % (cid, mol)) 
            for u in molatoms[mol]:
                mol1, atom1 = self.amp.get_mol_and_atom(u)
                at1 = self.amp.atom_types[mol1][atom1]                
                if u in acolors:
                    color = acolors[u]
                else:
                    color = "black"
                of.write("    %s [label=\"%s%s\",shape=\"%s\",color=\"%s\"];\n" % (u, at1, u, atomshape, color))
            of.write("  };\n")
            cid += 1

        for r in reactions:
            re = self.amp.reactions[r]
            if draw_edge_labels:
                edge_label = r
            else:
                edge_label = ""
            if r in rcolors:
                color = rcolors[r]
            else:
                color = "black"
            for u in re.maps:
                if u not in atoms:
                    continue
                for v in re.maps[u]:
                    if v not in atoms:
                        continue
                    if u in atom_edge_colors and v in atom_edge_colors[u]:
                        color = atom_edge_colors[u][v]
                        
                    of.write("   %s->%s [color=\"%s\",label=\"%s\"]\n" % (u, v, color, edge_label))
        of.write("}\n")
        of.flush()
        self.call_dot(imagetype, of.name, "%s.%s" % (ofn, imagetype))
        of.close()

def test():
    amp = atommap.AtomMappingParser()
    amp.parse_reaction_dir(sys.argv[1])    
    agv = AtomGraphVisualizer(amp)

    res = ["R00011_1", "R00012_1", "R00013_1", "R00014_1", "R00015_1"]
    agv.draw_atom_graph(res, sys.argv[2])

if __name__ == "__main__":
    test()
            
