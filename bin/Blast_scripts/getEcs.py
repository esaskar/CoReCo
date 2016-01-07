#!/usr/bin/env python   # #!means interpreter can know this is executible, path: means the interpreter you want to use for python

import sys, re           # import package

reEC = re.compile("DE\s+EC=(.+);") # regulation expression patterns for reEC; \s means whitespace character;+ means more one or more replications; ( ) means sub ER and . means any character.

f = open(sys.argv[1])      #input uniprot.dat
o = open(sys.argv[2],"w")  #output file
id = ac = None
ecs = []
for s in f:
    if s.startswith("ID"):
        id = s.strip().split()[1] # .strip can delet the heading and trailing of str with blank; split can change string into list
    elif s.startswith("AC"):
        if ac == None:
            ac = s.split()[1].strip().strip(";")
        ## Else will overwrite the first row of ACs
    elif s.startswith("DE"):
        m = reEC.match(s)
        if m:
            ec1 = m.group(1)
            ec2 = ec1.split("{")[0].strip()
            ecs.append(ec2)
    elif s.startswith("//"):
	if len(ecs) == 0:
        	e = "?"
	else:
    		e = ",".join(ecs)
        print "%s\t%s\t%s" % (ac, id, e)
        o.write("%s\t%s\t%s\n" % (ac, id, e))
        id = ac = None
        ecs = []
