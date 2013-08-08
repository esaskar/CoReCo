Computes atom mappings using A* algorithm. Finds such a mapping between the reactant
and product side atoms, that the atom incurs minimal number of changes to the bonds
during the reaction. 

Works on standard KEGG mol files (no extensions) with reaction definitions as in 
file 'kegg_reactions.txt'. The reactions should be balanced and corresponding mol-files 
valid.

usage: java Mapper2000 -h

Requires lots of memory, use at least 4 gigs of heap size for kegg reactions.

If used or improved, please cite:

 Heinonen, M. and Lappalainen, S. and Mielikäinen, T. and Rousu, J. Computing Atom 
 Mappings for Biochemical Reactions without Subgraph Isomoprhism. J. Comp. Biol. 2011, 18:43-58




