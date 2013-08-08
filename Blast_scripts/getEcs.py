#!/usr/bin/env python   # #!means interpreter can know this is  		           executible, path: means the interpreter you 				   want to use for python

import sys, re           # import package

reEC = re.compile("DE\s+EC=(.+);") # regulation expression patterns   					     for reEC; \s means whitespace   					     character;+ means more one or 					     more replications; ( ) means sub 					     ER and . means any character.

f = open(sys.argv[1])      #input uniprot.dat
o = open(sys.argv[2],"w")  #output file
id = ac = None
ecs = []
for s in f:
    if s.startswith("ID"):
        id = s.strip().split()[1] # .strip can delet the heading and 						trailing of str with blank; 						split can change string into 						list
    elif s.startswith("AC"):
        ac = s.split()[1].strip().strip(";")
    elif s.startswith("DE"):
        m = reEC.match(s)
        if m:
            ecs.append(m.group(1))   #list can .append but str not; 						group(0) means the whole 						pattern find by reEC; group(1) 						means only add the pattern in 						() (here is (.+)).
    elif s.startswith("//"):
	if len(ecs) == 0:
        	e = "?"
	else:
    		e = ",".join(ecs)
        print "%s\t%s\t%s" % (ac, id, e)
        o.write("%s\t%s\t%s\n" % (ac, id, e))
        id = None
        ecs = []
