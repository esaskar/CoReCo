#!/usr/bin/env python

# Original script probably from Esa Pitkanen
# Merja Oja edited it to also read from stdin
# and to take into account that sometimes
# there are several lines of Accessions

#import sys, re
import sys, re, fileinput

reEC = re.compile("DE\s+EC=(.+);")

#f = open(sys.argv[1])

id = ac = None
ecs = []
#for s in f:
for s in fileinput.input():
    if s.startswith("ID"):
        id = s.strip().split()[1]
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
            #ecs.append(m.group(1))
    elif s.startswith("//"):
        print "%s\t%s\t%s" % (ac, id, ",".join(ecs))
        id = ac = None
        ecs = []
