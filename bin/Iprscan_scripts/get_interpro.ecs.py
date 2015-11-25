#!/usr/bin/env python

import sys,re
from map_ipr2ec import mapipr2ec

ec2go = open(sys.argv[1])  #open file--->ec2go.txt
ipr2go = open(sys.argv[2]) #open file--->ipr2go.txt
ipr2ec = open(sys.argv[3],"w") # save ipr2ec result

mapipr2ec(ec2go,ipr2go,ipr2ec)
