#!/usr/bin/env python
"""
Parse KEGG compound formulae from 'compound' file.
"""

import sys, datetime, re

f = open("%s/compound" % (sys.argv[1])) # kegg dir
o = open(sys.argv[2], "w")  # output

#reparens = re.compile("([A-Z][a-z]*\d*)+(\([A-Z][a-z]*\d*\)\d*[nxyz]*\-?\d?)+")
reformula = re.compile("([A-Z][a-z]*)(\d*)")

def parse_formula(s):
    m = reformula.findall(s)
    return ",".join(map(lambda x: "%s:%s" % (x[0], x[1] or "1"), m))

def valid(form):
    return not ("(" in form or ")" in form)

o.write("#Output of \"%s\" on %s\n" % (" ".join(sys.argv), datetime.datetime.now()))
for s in f:
    if s.startswith("ENTRY"):
        mol = s.strip().split()[1]
    elif s.startswith("FORMULA"):
        form = s.strip().split()[1]
        if not valid(form):
            o.write("%s\t?\t%s\n" % (mol, form))
            continue
        formula = parse_formula(form)
        o.write("%s\t%s\t%s\n" % (mol, formula, form))
