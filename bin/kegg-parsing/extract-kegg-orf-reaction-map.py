#!/usr/bin/env python

import sys

import common

kdir = sys.argv[1]  # kegg dir
ofn = sys.argv[2]   # output file
o = open(ofn, "w")

ec2r, r2ec = common.read_ec_list(open(common.FILE_EC_MAP))

target_species = sys.argv[3].split(",")  # comma-separated list of species or "-" for all

f = open("%s/enzyme" % (kdir))
geneson = 0
ec2gene = {}
for s in f:
    if s[0] != " ":
        geneson = 0
    if geneson == 0:
        if s.startswith("ENTRY"):
            ec = s.strip().split()[2]
            ec2gene[ec] = {}
        elif s.startswith("GENES"):
            geneson = 1
            species = s[12:15]
            if target_species == "-" or species in target_species:
                ec2gene[ec][species] = set(s[16:].strip().split())
        elif s.startswith("///"):
            #print ec
            ss = ""
            for sp in ec2gene[ec]:
                ss += "%s:%s;" % (sp, ",".join(list(ec2gene[ec][sp])))
            
            if ec in ec2r:
                rr = ",".join(ec2r[ec])
            else:
                rr = "?"
            o.write("%s\t%s\t%s\n" % (ec, rr, ss))
    else:
        ns = s[12:15].strip()
        if ns != "":
            species = ns
        if target_species == "-" or species in target_species:
            if species not in ec2gene[ec]:
                ec2gene[ec][species] = set()
            ec2gene[ec][species].update(s[16:].strip().split())

