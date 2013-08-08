#!/usr/bin/env python

import sys, os, shutil

def convert_enzyme(ddir, tdir):
    f = open("%s/enzyme" % (ddir))
    o = open("%s/enzyme" % (tdir), "w")
    call = cacc = 0
    sys.stdout.write("enzyme: ")
    for s in f:
        if s.startswith("ENTRY"):
            ec = s.split()[2]
            accept = "-" not in ec
            call += 1
            if accept: cacc += 1 
        if accept:
            o.write(s)
    print "%d/%d complete enzymes" % (cacc, call)

def convert_reaction(ddir, tdir):
    sys.stdout.write("reaction: ")
    f = open("%s/reaction" % (ddir))
    o = open("%s/reaction" % (tdir), "w")
    call = cacc = 0
    buf = ""
    gens = []
    for s in f:
        if s.startswith("ENTRY"):
            buf = ""
            rid = s[12:18]
            accept = True
        buf += s

        if "general" in s:
            accept = False

        if s.startswith("///"):
            if accept:
                o.write(buf)
                buf = None
                cacc += 1
            else:
                gens.append(rid)
            call += 1
    print "%d/%d non-general reactions" % (cacc, call)
    o.close()
    
    o2 = open("%s/general-reactions" % (tdir), "w")
    o2.write("\n".join(gens))
    o2.close()

def convert_compound(ddir, tdir):
    sys.stdout.write("compound: ")
    shutil.copy("%s/compound" % (ddir), "%s/compound" % (tdir))
    print "copied"

def main(ddir, tdir):
    convert_enzyme(ddir, tdir)
    convert_reaction(ddir, tdir)
    convert_compound(ddir, tdir)

if __name__ == "__main__":
    ddir = sys.argv[1]  # kegg directory
    tdir = sys.argv[2]  # target directory
    main(ddir, tdir)
