#!/usr/bin/env python   # 

import sys,fileinput,re,os         # import package

reExtracomma=re.compile("'")

#print(sys.argv)
print "Combining EC numbers"

if len(sys.argv) >= 9:     			
    fECsprot = open(sys.argv[1])        # uniprot_sprot.ec_files.txt
    fEComat =  open(sys.argv[2]) 	# gene2eclist.txt
    if os.path.exists( sys.argv[3] ):
    	fidMap_add =  open(sys.argv[3]) # additional_Uniprot_ids.txt		
    if os.path.exists( sys.argv[4] ):
        fidMap =  open(sys.argv[4])     # gene2eclist_Uniprot_ids.txt
    ## If both of above missing then assuming that all genes in gene2eclist already have uniprot IDs
    ## also works in sys.arg[4] is a file of size 0
    output =  open(sys.argv[5],"w") 	# output -> ec_files.txt
    output2 = open(sys.argv[6],"w") 	# output -> not_in_swissprot.ids    
    output3 = open(sys.argv[7],"w")	# overlapping table	
    output4 = open(sys.argv[8],"w") 	# temp	
    if len(sys.argv) > 9:
        appendornot =  sys.argv[9] 	        # if sys.arv=="append", append ECs instead of replace by ECs from models (default)
        print("Appending ECs")
 	
#    EComat =  fileinput.input()
else:
    print "Not enough inputs"

output3.write("%s\t%s\t%s\t%s\n" % ("SwissProt gene ID", "SwissProt EC" , "Users gene ID" , "Users EC"))

id = ac = None
ecs = []
missing = []

oma2uniprot = {}
uniprot2oma = {}

# Adds the additional Uniprot Ids that can't be found with the JavaTools program. These additional Uniprot Ids should be put in a tab separated file called "additional_Uniprot_ids.txt"
if os.path.exists( sys.argv[3] ):
	for s in fidMap_add:
		apu = s.strip("\n").split("\t")
		if len(apu) == 2:	
			id, acc = apu	
			if acc != "null":
				oma2uniprot[id]=acc
				uniprot2oma[acc]=id
			else: 
				missing.append(apu)

# Adds the Uniprot Ids found with the JavaTools program
if os.path.exists( sys.argv[4] ):
    for s in fidMap:
        apu = s.strip("\n").split(", ")
        if len(apu)==2:
            id, acc = apu
            if reExtracomma.match(apu[1]):
                acc = apu[1].split("'").split("_");
                print id
                print acc
            if acc != "null":
                oma2uniprot[id]=acc
                uniprot2oma[acc]=id
            else: 
                missing.append(apu)
        elif len(apu)==3:
            id, acc, puppu = apu
            if reExtracomma.match(apu[1]):
                acc = apu[1].split("'").split("_");
                print id
                print acc
            if acc != "null":
                oma2uniprot[id]=acc
                uniprot2oma[acc]=id
            else: 
                missing.append(apu)
        else:        
            missing.append(apu)

#print(missing)

uniprotEC = {}
for s in fECsprot:
    id,acc,ecs = s.strip("\n").split("\t")
    uniprotEC[id] = [acc, ecs]

notinUniprotEC = {}
omatECfound = []
for s in fEComat:
    id,acc,ecs = s.strip("\n").strip("\r").split("\t")
    id = id.split(".")[0]
    
    if id in oma2uniprot:
        idUniprot = oma2uniprot[id]
    else:
        if id not in missing:
            print "id " + id + " not recognized, assuming it is an uniprotID"
        idUniprot = id

    if idUniprot!="":
        if idUniprot in uniprotEC:
            accUniprot, ecsUniprot = uniprotEC[idUniprot]
            output3.write("%s\t%s\t%s\t%s\n" % (accUniprot, ecsUniprot, acc, ecs))
            #print("%s\t%s\t%s\t%s" % (accUniprot, ecsUniprot, acc, ecs))
            if appendornot=="append":
                eclistOma = list(ecs.split(","))
                eclistUniprot = list(ecsUniprot.split(","))
                for ecoma in eclistOma: 
                    if ecoma not in eclistUniprot:
                        eclistUniprot.append(ecoma)
                if "" in eclistUniprot:
                    print(eclistUniprot)
                    eclistUniprot.remove("")
                    print(eclistUniprot)
                #print(eclistUniprot)
                ecs = ",".join(eclistUniprot)
            #else:
            #    print("\n")
            #    ## Brutally replace information that was in Uniprot with the info from the model
            #    ## would be possible to split both and combine and keep unique
            uniprotEC[idUniprot] = [accUniprot,ecs]
        else:
            ## Add new field to Uniprot
            uniprotEC[idUniprot] = [acc,ecs]
            ## Add identifier to the list for which sequence data will be retrieved
            omatECfound.append(idUniprot)
    else:
        print("idUniprot was empty: %s" % (idUniprot))
        notinUniprotEC[id] = [acc,ecs]
    


for id in uniprotEC:
    output.write("%s\t%s\t%s\n" % (id,uniprotEC[id][0],uniprotEC[id][1]))

for id in omatECfound:
    output2.write("%s\n" % (id))

for id in notinUniprotEC:
    output4.write("%s\t%s\t%s\n" % (id,notinUniprotEC[id][0],notinUniprotEC[id][1]))




