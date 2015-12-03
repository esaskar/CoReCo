#!/usr/bin/python
# -*- coding: iso-8859-15 -*-


import string, glob, optparse, sys, os, re
import kegg
from MPfunctions import *
from MessagePassing import *


ATOM_NUMBER = 0
FEATURE_NAME = 1
CONTEXT = 2
VALUE = 3

# features are encoded as strings
#   MOL:ATOM_NUMBER:FEATURE_NAME:CONTEXT=VALUE
# features = { mol:[atom_nmbr, feature_name, context, value], mol:[...],..}
features = {}
molfiles = []
mols = []
context = []


#############
#
#  args
#

parser = optparse.OptionParser(version = "1.00", description = "Generates local features", usage = "%prog [-something] MOLFILES")

parser.add_option("-v", "--verbose", dest="verbose", action="store_true", help="Verbose mode")
parser.add_option("-q", "--quiet", dest="quiet", action="store_true", default=False, help="Quiet mode")
parser.add_option("-d", "--debug", dest="debug", action="store_true", default=False, help="Debug mode")
parser.add_option("-i", "--inputfile", dest="inputfile", default="", help = "list of mols")
parser.add_option("-m", "--moldir", dest="moldir", default="", help = "directory where mol files are located")
parser.add_option("-o", "--output-dir", dest="output", default="", help="Output folder")
parser.add_option("-k", "--context", dest="context", type="string", default="0-1000", help="Context size, <int> or <lb>-<ub>, e.g. '-k 3' or '-k 0-2'. Defaults to all possible contexts.")
parser.add_option("-x", "--exclude", dest="exclude", action="store_true", default=False, help="Exclude central atom info from feature")
parser.add_option("-a", "--full-atominfo", dest="atominfo", action="store_true", default=False, help="Use full atominfo")
parser.add_option("-b", "--full-bondinfo", dest="bondinfo", action="store_true", default=False, help="Use full bondinfo")
parser.add_option("-n", "--noappend", dest="noappend", action="store_true", default=False, help="Use full bondinfo")

featuregroup = optparse.OptionGroup(parser, "Features")
featuregroup.add_option("-A", dest="ad", action="store_true", help="Atom distribution")
featuregroup.add_option("-B", dest="bd", action="store_true", help="Bond distribution")
featuregroup.add_option("-W", dest="wiener", action="store_true", help="Wiener index")
featuregroup.add_option("-R", dest="ring", action="store_true", help="Ring")
featuregroup.add_option("-M", dest="morgan", action="store_true", help="Morgan index")
parser.add_option_group(featuregroup)


(options, args) = parser.parse_args()


if not args and not options.inputfile:
	parser.error("No mol files defined")

if "-" in options.context:
	s = options.context
	context = range(int(s.split("-")[0]), int(s.split("-")[1])+1)
else:
	context = [int(options.context.strip())]

if not context:
	parser.error("No context given")

mols = args

if options.inputfile:
	f = open(options.inputfile)
	for line in f:
		line = line.strip()
#		if options.moldir:
#			molfiless.append(options.moldir + line + ".mol")
#		else:
#			molfiless.append(line + ".mol")
		mols.append(line)
	f.close()

if options.quiet:
	sys.stdout = open(os.devnull, "w")


def WriteResults(mols, features):
	# convert molfile into featfile
	for mol in mols:
		
		if options.verbose:
			print "writing results for ", mol
		
		molfile = ""

		if options.moldir:
			molfile = options.moldir + mol + ".mol"
		else:
			molfile = mol + ".mol"
		
		featfile = ""
		
		if options.output:
			featfile = options.output + mol + ".feat"
		else:
			featfile = mol + ".feat"
		
		if options.verbose:
			print "Mol-file", molfile, "with feature-file", featfile
		
		existing_lines = []
		
		# file exists
		if os.path.isfile(featfile):
			ff = open(featfile)
			existing_lines = ff.read()
			ff.close()
			
			if options.verbose:
				print "Feat-file already exists"
		
		f = open(featfile, "a")
		
		for feat in features[mol]:
			s = "%s:%d:%s:%d=%s" % (mol, feat[ATOM_NUMBER], feat[FEATURE_NAME], feat[CONTEXT], feat[VALUE])
			if s not in existing_lines:
				if options.verbose:
					print "writing", s
				f.write(s + "\n")
			else:
				if options.verbose:
					print "line", s, "exists, skipping"
		
		f.close()
		
		print "Mol", mol, "written to", featfile


def ExtractValues(graph):
	data = {}
	for n in graph.atoms:
		data[n.id] = n.val[-1]
	return data


def ExtractFeatures(mol, graph, feature_name, context, join=True, sort=True):
	if not features.has_key(mol):
		features[mol] = []
	
	s = ""
	for n in graph.atoms:
		contextset = [i for i,v in enumerate(n.val) if i in context]

		for c in contextset:
			data = None
			try:
				data = n.val[c].vec
				if sort:
					data.sort()
			except:
				data = n.val[c]
			
			if join:
				features[mol].append( [n.id, feature_name, c, ''.join(data)] )
			else:
				features[mol].append( [n.id, feature_name, c, data] )
			
		
		if not contextset and feature_name == "RING":
			for c in context:
				features[mol].append( [n.id, feature_name, c, "False"] )


def ComputeFeatures():
	mfs = len(mols)
	print mfs, "molfiles ..."
	
	for mol in mols:
		
		print "computing features for ", mol
		
		molfile = ""

		if options.moldir:
			molfile = options.moldir + mol + ".mol"
		else:
			molfile = mol + ".mol"

		if options.output:
			featfile = options.output + mol + ".feat"
		else:
			featfile = mol + ".feat"

		if os.path.isfile(featfile) and options.noappend:
			print "featurefile already exists for ", mol
			print "not recomputing "
			continue
		
#		if options.verbose:
		print "Molfile", molfile

		if os.path.exists(molfile):
			m = kegg.Molecule(mol, molfile=molfile)
			g = m.graph
			
			# if atom symbols contain extra information, strip it if not wanted
			if not options.atominfo:
				m.strip_atominfo()
				
				if options.bondinfo:
					m.full_bondinfo()
		
			# for each kegg mol-file, call the algorithms by feature wanted
			# call the message-passing algorithm with appropriate sumfunctions as parameters
			# e.g. for atom distribution, give as parameters the functions:
			#  - atom_node_feat function to simply give the atoms label as its basic value
			#  - msg_feat function to simply set the value of msg the same as behind node's
			#  - atom_distribution_node_sumfunc function to sum up all incoming messages 
			#    (containing labels)
			#  - atom_distribution_msg_sumfunc function to set msgs as "everything-behind"
		
			# The algorithmB incovations put the results directly into the nodes (Graph.Node), 
			# retrieve the wanted features using ExtractFeatures() function
		
			if options.verbose:
				print "computing feats for", molfile
		
			# WARNING: doesn't work correctly, the algorithmB needs different form for edge based message passing...
			if options.bd:
				AlgorithmB(g, bond_node_feat, bond_msg_feat, bond_distribution_node_sumfunc, bond_distribution_msg_sumfunc, verbose=options.debug, maxiter=max(context))
				ExtractFeatures(mol, g, "BOND_DISTRIBUTION", context, True)
		
			if options.ad and not options.exclude:
				AlgorithmB(g, atom_node_feat, msg_feat, atom_distribution_node_sumfunc, atom_distribution_msg_sumfunc, verbose=options.debug, maxiter=max(context))
				ExtractFeatures(mol, g, "ATOM_DISTRIBUTION", context, True)
		
			if options.ad and options.exclude:
				AlgorithmB(g, atom_node_feat, msg_feat, atom_distribution_exclude_node_sumfunc, atom_distribution_msg_sumfunc, verbose=options.debug, maxiter=max(context))
				ExtractFeatures(mol, g, "ATOM_DISTRIBUTION_EXCLUDE", context, True)
		
			if options.ring:
				AlgorithmB(g, distance_node_feat, msg_feat, ring_node_sumfunc, min_distance_msg_sumfunc, verbose=options.debug, maxiter=max(context))
				ExtractFeatures(mol, g, "RING", context, False)
		
			if options.wiener:
				AlgorithmB(g, DistanceNodeFeat(), msg_feat, min_distance_node_sumfunc, min_distance_msg_sumfunc, verbose=options.debug, limit_rounds=True)
			
				# fs is a dict
				fs = ExtractValues(g)
				AlgorithmB(g, ExtdataNodeFeat(fs), msg_feat, wiener_node_sumfunc, atom_distribution_msg_sumfunc, verbose=options.debug, maxiter=max(context))
			
				ExtractFeatures(mol, g, "WIENER", context, False, False)
		
			if options.morgan:
				AlgorithmB(g, MorganNodeFeat(), MorganMsgFeat(), MorganNodeSumfunc(), MorganMsgSumfunc(), verbose=options.debug, maxiter=max(context), selfloops=True)
				ExtractFeatures(mol, g, "MORGAN", context, False, False)
		
		
			print "Mol", mol, "computed"

			WriteResults([mol], features)

		else:
			sys.stderr.write("Cannot find molfile " + molfile + "\n")
	


if __name__ == "__main__":
	ComputeFeatures()
#	WriteResults(molfiles, features)




