#!/usr/bin/python
# -*- coding: iso-8859-15 -*-
#
#
# parsers for:
# - kegg ligand reaction file
# - kegg ligand compound file
# - kegg ligand mol files
#


def parse_reaction(filename="/home/group/icomic/data/kegg/ligand/LATEST/reaction"):
	f = open(filename)
	lines = f.readlines()
	f.close()
	
	D = {}
	r = None  # current reaction id
	block = None  # current block id
	blockdata = []
	
	for i,line in enumerate(lines):
		prefix = line[0:12].strip()
		data = line[12:].rstrip()
		
		if prefix:
			if block:
				if block == "ENTRY":
					r = blockdata[0].split()[0]
					D[r] = {}
#					print r
				elif block == "NAME":
					D[r][block] = " ".join(blockdata)
				elif block == "DEFINITION":
					D[r][block] = " ".join(blockdata)
				elif block == "EQUATION":
					D[r][block] = " ".join(blockdata)
				elif block == "RPAIR":
					D[r][block] = blockdata
				elif block == "ENZYME":
					D[r][block] = [s for row in blockdata for s in row.split()]
				elif block == "PATHWAY":
					D[r][block] = blockdata
				elif block == "ORTHOLOGY":
					D[r][block] = " ".join(blockdata)
				elif block == "COMMENT":
					D[r][block] = " ".join(blockdata)
				elif block == "REFERENCE":
					D[r][block] = " ".join(blockdata)
				elif block == "REMARK":
					D[r][block] = " ".join(blockdata)
				elif block == "///":
					pass
				else:
					print "unknown block:", block
			
			block = prefix
			blockdata = [data.strip()]
		else:
			if data.startswith(" "):
				blockdata[-1] += " " + data.strip()
			else:
				blockdata.append(data)
	
	return D



#R = parse_reaction("/home/mqheinon/temp/reaction")


# okei... eli: puretaan uusin ligand icomiciin, muodostetaan lista valideista reaktioista käyttäen
# pythonin Reaction ja Molecule luokkia...
# pidetään vedyt poissa kokonaan
# sillä selvä!








