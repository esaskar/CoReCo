#!/usr/bin/evn python

import sys,re


def mapipr2ec(ec2go,ipr2go,ipr2ec):
        #save ec2go in gos, save ipr2go in iprs and save ec in ecs
        iprs ={}
        ecs = []
        gos = {} 

        for s in ipr2go:  # extract all ipr2go and save it into a dictionary iprs; iprs[ipr] = [ipr, GO]
                if len(s.strip().split("\t"))>1: # exclude those iprs without GO term.
                        ipr, GO = s.strip().split("\t")
                #       print ipr,GO
                        iprs[ipr] = [ipr, GO]
                else:
                        ipr = s.strip().split("\t")[0]  # save those ipr without GO term.
                #       print ipr
                        iprs[ipr] = [ipr, '?']
                #print iprs[ipr][1]

        for f in ec2go:
                if f.startswith("EC"):    #EC:1.1.1 > GO:oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor ; GO:0016616
                        ec = f.strip().split(">")[0].strip().split(":")[1].strip()  #extract EC files :   1.1.1
                        go_id = f.strip().split(";")[1].strip()
                #       print ec, go_id  # ec and go_id looks like EC:6.3.4.7 GO:0050260
                        gos[go_id] = [ec]
                #       print gos[go_id]
        gokeys = gos.keys()
        gokeys.sort()
        iprkeys = iprs.keys()
        iprkeys.sort()

        for ipr in iprkeys:
                GO = iprs[ipr][1]  # GO looks like GO:0000724,GO:0019789,GO:0030915,GO:0034184,GO:0045842 or ?
                # split GO into list that containing a list of go
                GOs = GO.strip().split(",")[:]
                for go in GOs:
                        for gokey in gokeys:
                                if go==gokey:
                                        ecs.append(gos[gokey][0])
                        #               print ecs
        #       print ecs
                if len(ecs) == 0:
                        ec = '?'
                else:
                        ec = ",".join(ecs)
        #       print ec
        #       print "%s\t%s\t%s" % (iprs[ipr][0], iprs[ipr][1], ec)
                ipr2ec.write("%s\t%s\t%s\n" % (iprs[ipr][0], iprs[ipr][1], ec))
                ecs = []

if __name__ == "__main__":
        #map ipr2ecs by ipr2go and ec2go
        ec2go = open(sys.argv[1])  #open file--->ec2go.txt
        ipr2go = open(sys.argv[2]) #open file--->ipr2go.txt
        ipr2ec = open(sys.argv[3],"w") # save ipr2ec result
        mapipr2ec(ec2go,ipr2go,ipr2ec)

