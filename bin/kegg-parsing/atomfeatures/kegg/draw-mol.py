#!/usr/bin/python
#




import kegg
import optparse, sys




parser = optparse.OptionParser(usage="%prog [opts] molecule")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Verbose mode")
parser.add_option("-q", "--quiet", dest="quiet", action="store_true", default=False, help="Quiet mode")

(options, args) = parser.parse_args()

if not args:
	parser.error("no molecule given")

mol = kegg.Molecule(args[0])

md = kegg.MolDraw(mol)
md.layout()
mol.svg()

print "done"









