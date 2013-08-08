#!/usr/bin/env python

import sys, tempfile, os, subprocess

import common, plot

def splot(datafn, plotfn, format):
    tf2 = tempfile.NamedTemporaryFile()
    o2 = open(tf2.name, "w")

    o2.write(plot.plot_format("%s.scatter" % (plotfn), format))

    o2.write("""
d = as.matrix(read.table("%s"))
par(mfrow=c(3,1), oma=c(2, 0, 3, 0), omi=c(0, 0, 0.8, 0))
plot(d[,2], d[,1], ylab="Posterior probability", xlab="BLAST score")
plot(d[,3], d[,1], ylab="Posterior probability", xlab="GTG score")
plot(d[,2], d[,3], xlab="BLAST score", ylab="GTG score")
dev.off()
""" % (datafn))

    o2.write(plot.plot_format("%s.density" % (plotfn), format))
    o2.write("""
par(mfrow=c(3,1), oma=c(2, 0, 3, 0), omi=c(0, 0, 0.8, 0))
plot(density(d[,1], from=0, to=1), xlab = "Posterior probability", main = "")
plot(density(d[,2], from=0, to=max(d[,2])), xlab = "BLAST score", main = "")
plot(density(d[,3], from=0, to=max(d[,3])), xlab = "GTG score", main = "")
""")

    o2.write(plot.plot_format("%s.full-vs-naive" % (plotfn), format))
    o2.write("""
par(mfrow=c(2,1), oma=c(2, 0, 3, 0), omi=c(0, 0, 0.8, 0))
plot(d[,1], d[,4], xlab = "Posterior, full model", ylab = "Posterior, naive model", main = "")
plot(density(d[,5], from=min(d[,5]), to=max(d[,5])), xlab = "Log ratio, full vs naive", main = "")
""")

    o2.flush()
    subprocess.call("R CMD BATCH %s" % (tf2.name), shell = True)
    o2.close()

def plot_reaction_scores(ddir):
    scores, all_ecs = common.read_scores(ddir, augmented_format = True)

    keys = scores.keys()
    keys.sort()
    try:
        os.mkdir("%s/%s" % (ddir, common.REACTION_SCORE_PLOT_DIR))
    except:
        pass

    for species in keys:
        sys.stdout.write("%s " % (species))
        plotfn = "%s/%s/%s" %(ddir, common.REACTION_SCORE_PLOT_DIR, species)
        tf = tempfile.NamedTemporaryFile()
        o = open(tf.name, "w")
        for ec in scores[species]:
            rscore = scores[species][ec]
            o.write("%s\n" % (rscore))
        o.flush()
        splot(tf.name, plotfn, "pdf")
        splot(tf.name, plotfn, "png")
        o.close()

if __name__ == "__main__":
    ddir = sys.argv[1]
    plot_reaction_scores(ddir)
