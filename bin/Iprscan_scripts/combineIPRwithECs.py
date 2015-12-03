#!/usr/bin/env python

import sys,re


header = ["Org", "SeqId", "IPRScanMethod", "MethodRef", "EValue", "iprID", "GO", "ECs"]


ecf = open(sys.argv[1]) # interpro.ecs from ipr2ec.py
rawf = open(sys.argv[2]) # interpro raw result file
sof = open(sys.argv[3])  # seq_org_list
result = open(sys.argv[4],"w") # save combineIPRwithECs results
iprs = {}


#wirte title in the result
result.write("%s\n" % "\t".join(header))
# save each raw in a dictionary like: iprs[iprid] = [ipr, GO, ecs]
for f in ecf:
	iprid, GO, ecs = f.strip().split("\t")
	iprs[iprid] = [iprid, GO, ecs]
#['IPR026873', 'GO:0004663,GO:0018344', '2.5.1.60']
#	print iprs[iprid]

seqtoorg = {}
for s in sof:
	seqid, orgid = s.strip().split("\t")
	seqtoorg[seqid] = orgid
#seem:	Psti PAS_FragD_0013
#	print seqtoorg[seqid], seqid


for s in rawf:
	vals = s.strip().split("\t")
#	print len(vals), vals
	seqid, checksum, seqlen, method, mres1, mres2, matchstart, matchend, evalue, matchstatus, date, iprid, iprdesc1 = vals[0:13]
	if len(vals) == 14:
		iprdesc2 = vals[13]
	else:
		iprdesc2 = None

	if seqid not in seqtoorg:
		sys.stderr.write("Cannot find seq %s in org list\n" % (seqid))
		org = "?"
	else:
		org = seqtoorg[seqid]
    #there is something wrong
	if iprid != "NULL" and iprid in iprs:
		ip = "\t".join(iprs[iprid])
		
    	else:
		ip = "%s%s" % (iprid, "\t?" * (len(["iprID", "GO", "ECs"]) - 1))
	result.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (org, seqid, method, mres1, evalue, ip))


