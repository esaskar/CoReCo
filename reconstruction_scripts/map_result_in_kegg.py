#!/usr/bin/env python

import sys, shutil


sys.path.append("../model_training_scripts/")
import common

KEGG_PATH_LINK = "http://www.genome.jp/kegg-bin/show_pathway?%s"

JS_SORTABLE_TABLE = "webtoolkit.sortabletable.js"
HTML_HEAD_SORTABLE_TABLE = """
  <script type="text/javascript" src="%s"></script>
        <style>
                table {
                        text-align: left;
                        font-size: 12px;
                        font-family: verdana;
                        background: #c0c0c0;
                }
 
                table thead  {
                        cursor: pointer;
                }
 
                table thead tr,
                table tfoot tr {
                        background: #c0c0c0;
                }
 
                table tbody tr {
                        background: #f0f0f0;
                }
 
                td, th {
                        border: 1px solid white;
                }
        </style>
""" % (JS_SORTABLE_TABLE)

HTML_BODY_SORTABLE_TABLE = """
  <script type="text/javascript">
  var t = new SortableTable(document.getElementById('%s'), 100);
  </script>
"""

def main(fillsfn, odir):
    fills = common.read_fills(fillsfn) # reco res dir

    #f = open("%s/fills" % (sys.argv[1]))  # reconstruction result dir
#    f2 = open("aux/reaction-pathways")  # reaction -> pathway1,pathway2,...
#    f3 = open("aux/pathway-names") # pathway id -> pathway name

    f2 = open("../../data/Kegg/aux/reaction-pathways")  # reaction -> pathway1,pathway2,...
    f3 = open("../../data/Kegg/aux/pathway-names") # pathway id -> pathway name

    R = {}
    all_fills = set()

    for r in fills:
        for r2 in fills[r].reactions:
            base = r2.split("_")[0]
            if base != r:
                all_fills.add(base)

    path_r = {}
    for s in f2:
        vals = s.strip().split()
        paths = []
        if len(vals) > 1:
            paths = vals[1].split(",")
        for p in paths:
            if p not in path_r:
                path_r[p] = set()
            path_r[p].add(vals[0])

    path_names = {}
    for s in f3:
        pid, name = s.strip().split("\t")
        path_names[pid] = name

    o = open("%s/index.html" % (odir), "w")
    o.write("""
    <html>
    <head>
    <title>Pathways from fills</title>
    %s
    </head>
    <body>
    <p>
    Reconstruction result in %s.
    </p>
    <table id=\"restable\">
    <thead>
    <tr>
    <th class="c1">Pathway</th>
    <th class="c2">Reactions</th>
    <th class="c3">NotGapped</th>
    <th class="c4">GapFilled</th>
    <th class="c5">Filler</th>
    <th class="c6">FillTooCostly</th>
    <th class="c7">NoGapFixFound</th>
    <th class="c8">Unknown</th>
    <th class="c9">NotInReconstruction</th></tr>
    </thead>
    """ % (HTML_HEAD_SORTABLE_TABLE, sys.argv[1]))

    o.write("<tbody>\n")
    for path in path_r:
        pid = "rn%s" % (path[2:])
        link_s = KEGG_PATH_LINK % (pid)

        link_s += "/default%3d%23fefefe"

        if pid in path_names:
            pname = path_names[pid]
        else:
            pname = pid

        nr = len(path_r[path])
        notgapped = costly = filled = nofix = unknown = notinreco = filler = 0
        all_gaps = set()
        for r in path_r[path]:
            fgcolor = "000000"
            bgcolor = "fefefe"
            all_gaps.add(r)
            if r in fills:
                dec = fills[r].decision
                #all_fills.update(fills)
                if dec == "NotGapped":
                    bgcolor = "0000fe"
                    notgapped += 1
                elif dec == "FillTooCostly":
                    bgcolor = "fefe00"
                    costly += 1
                elif dec == "GapFilled":
                    bgcolor = "00fe00"
                    filled += 1
                elif dec == "NoGapFixFound":
                    bgcolor = "fe0000"
                    nofix += 1
                else:
                    bgcolor = "aaaaaa"
                    unknown += 1
            else:
                if r in all_fills:
                    filler += 1
                    bgcolor = "00fefe"
                else:
                    bgcolor = "fefefe"
                    #bgcolor = None
                    notinreco += 1
            #link_s += "/%s" % (r)
            if bgcolor != None:
                link_s += "/%s%%09%%23%s" % (r, bgcolor)
            #link_s += "/%s%%09%%23%s,%%23%s" % (r, bgcolor, fgcolor)

        #rest = all_fills.difference(all_gaps).intersection(path_r[path])
        #rest = all_fills.difference(all_gaps)
        rest = all_fills.intersection(path_r[path])

        #print rest

        bgcolor = "00fefe"
        for r in rest:
            link_s += "/%s%%09%%23%s" % (r, bgcolor)

        o.write("<tr><td class=\"c1\"><a href=\"%s\">%s</a></td><td class=\"c2\">%s</td><td class=\"c3\">%s</td><td class=\"c4\">%s</td><td class=\"c5\">%s</td><td class=\"c6\">%s</td><td class=\"c7\">%s</td><td class=\"c8\">%s</td><td class=\"c8\">%s</td></tr>\n" % (link_s, pname, nr, notgapped, filled, filler, costly, nofix, unknown, notinreco))

    o.write("""
    </tbody>
    </table>
    %s
    </body>
    </html>
    """ % (HTML_BODY_SORTABLE_TABLE % ("restable")))

    print "Wrote %s/index.html" % (odir)

    try:
        shutil.copy(JS_SORTABLE_TABLE, "%s/%s" % (odir, JS_SORTABLE_TABLE))
    except IOError:
        print "Warning: unable to copy %s to %s" % (JS_SORTABLE_TABLE, odir)

if __name__ == "__main__":
    fillsfn = sys.argv[1] # reco res dir
    odir = sys.argv[2]  # output dir
    main(fillsfn, odir)

