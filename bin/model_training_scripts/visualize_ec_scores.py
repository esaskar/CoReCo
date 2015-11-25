#!/usr/bin/env python

import sys, os, tempfile, subprocess

import common

#digraph structs { node [shape=record]; 
#struct1 [label="<f0> left|<f1> mid\ dle|<f2> right"]; 
#struct2 [label="<f0> one|<f1> two"]; 
#struct3 [label="hello\nworld |{ b |{c|<here> d|e}| f}| g | h"]; struct1:f1 -> struct2:f0; struct1:f2 -> struct3:here; } 

#struct2 [label=< <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0"> <TR><TD PORT="f0">one</TD><TD>two</TD></TR> </TABLE>>]; 

def name_to_label(x):
    return "node_%s" % (x)

def score_to_color(x):
    g = x
    r = 1.0 - g
    b = 0.25
    #print r, g, b, "#%.2x%.2x%.2x" % (r * 255, g * 255, b * 255)
    return "#%.2x%.2x%.2x" % (r * 255, g * 255, b * 255)

def norm_blast(x, max_blast):
    return x / max_blast

def tree_to_dot(tree, scores, max_blast, ec):
    parent, children = tree
    s = "digraph ectree {\nnode [shape=record];\n\tlabel=< <FONT POINT-SIZE=\"32\">EC %s</FONT> >;\nlabelloc=\"t\";\nlabelfontsize=24;\n" % (ec)

    for u in children:
        #print u

        leaf = len(children[u]) == 0
        sc = scores[u]

        fullp = sc.pscore
        naivep = sc.npscore
        blast = sc.bscore
        gtg = sc.gscore
        fullcol = score_to_color(fullp)
        naivecol = score_to_color(naivep)
        blastcol = score_to_color(norm_blast(blast, max_blast))
        gtgcol = score_to_color(gtg)
        if not leaf:
            s += "\t%s [label=< <TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\"><TR><TD COLSPAN=\"2\">%s</TD></TR><TR><TD COLSPAN=\"2\" BGCOLOR=\"%s\">%.3f</TD></TR></TABLE>>];\n" % (name_to_label(u), u, fullcol, fullp)
        else:
            s += "\t%s [label=< <TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\"><TR><TD COLSPAN=\"2\">%s</TD></TR><TR><TD COLSPAN=\"2\" BGCOLOR=\"%s\">%.3f</TD></TR><TR><TD COLSPAN=\"2\" BGCOLOR=\"%s\">%.3f</TD></TR><TR><TD COLSPAN=\"1\" BGCOLOR=\"%s\">%d</TD><TD COLSPAN=\"1\" BGCOLOR=\"%s\">%.2f</TD></TR></TABLE>>];\n" % (name_to_label(u), u, fullcol, fullp, naivecol, naivep, blastcol, blast, gtgcol, gtg)
    for u in children:
        for v in children[u]:
            s += "\t%s -> %s;\n" % (name_to_label(u), name_to_label(v))
    s += "}\n"
    return s

def plot(dots, outdir, outfn, format = "png"):
    of = tempfile.NamedTemporaryFile()
    of.write(dots)
    of.flush()
    subprocess.call("dot -T%s -o %s/%s.%s %s" % (format, outdir, outfn, format, of.name), shell = True)
    of.close()

def main(ddir):
    tree = common.read_tree_adj_list(ddir)
    scores, all_ecs = common.read_scores(ddir, augmented_format = True)

    ecscores = {}
    max_blast = 0
    for species in scores:
        for ec in scores[species]:
            if ec not in ecscores:
                ecscores[ec] = {}
            ecscores[ec][species] = scores[species][ec]
            blast = scores[species][ec].bscore
            if blast > max_blast:
                max_blast = blast

    try:
        os.mkdir("%s/%s" % (ddir, common.EC_VISUALS_DIR))
    except:
        pass

    keys = ecscores.keys()
    keys.sort()
    for ec in keys:
        #print ec
        s = tree_to_dot(tree, ecscores[ec], max_blast, ec)
        #print s
        outdir = "%s/%s" % (ddir, common.EC_VISUALS_DIR)
        outfn = "%s" % (ec)
        plot(s, outdir, outfn)
        #sys.exit(1)

if __name__ == "__main__":
    ddir = sys.argv[1] # project dir
    main(ddir)
