#!/usr/bin/env python

import sys

from metabolism.parser import kegg_ligand_parser as klp

f = open(sys.argv[1]) # kegg compound file
o = open(sys.argv[2], "w")  # output
voc = klp.parse_vocabulary_from_compound(f)

ids = list(voc.ids)
ids.sort()
for mid in ids:
    o.write("%s\t%s\t%s\n" % (mid, voc.getCommonName(mid), voc.getNonIdSynonym(mid)))

