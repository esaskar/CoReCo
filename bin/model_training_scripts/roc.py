#!/usr/bin/env python

import sys, os, tempfile, subprocess

import common, plot

FILTER_PARTIAL_ECS = 0

def filter_partial_ecs(ecs):
    A = set()
    for ec in ecs:
        if "-" not in ec:
            A.add(ec)
    return A
            
def compute_roc_curve(models, scores, all_ecs, of, species, index, max_value, filter_partials = False):
    # argument index gives the score attribute name in Score
    num_pos = len(models[species])
    num_neg = len(all_ecs) - num_pos

    tp = {}  # prob -> tp
    fp = {}  # prob -> fp

    ecs = scores[species].keys()

    # sort by posterior probability in full model or naive model
    ecs.sort(lambda x, y: -cmp(scores[species][x].__getattribute__(index), scores[species][y].__getattribute__(index)))

    if filter_partials:
        ecs = filter(lambda x: "-" not in x, ecs)

    pos = 0  # num positives, cumulative
    neg = 0  # num negatives, cumulative

    if max_value == None:
        max_value = scores[species][ecs[0]].__getattribute__(index)

    for i, ec in enumerate(ecs):
        p = scores[species][ec].__getattribute__(index)
        if ec in models[species]:
            pos += 1
        else:
            neg += 1
        tp[p] = pos
        fp[p] = neg
        #print i, p, pos, neg

    k = tp.keys()
    k.sort(lambda x, y: cmp(x, y))
    of.write("#Param FPR TPR F1 Prec AUC\n")
    of.write("0.0\t1.0\t1.0\t0.0\t0.0\t0.0\n")
    pfpr = 1.0
    ptpr = 1.0
    auc = 0.0
    for p in k:
        fn = num_pos - tp[p]
        prec = 1.0 * tp[p] / (tp[p] + fp[p])
        #sens = 1.0 * tp[p] / (tp[p] + fn)
        tpr = 1.0 * tp[p] / num_pos #TPR/sensitivity
        assert(tpr >= 0 and tpr <= 1)
        fpr = 1.0 * fp[p] / num_neg #FPR
        assert(fpr >= 0 and fpr <= 1)
        if prec + tpr == 0:
            f1 = 0
        else:
            f1 = 2 * prec * tpr / (prec + tpr)

        auc += (pfpr - fpr) * tpr + (pfpr - fpr) * (ptpr - tpr) / 2
        #ccc = (pfpr - fpr) * tpr + (pfpr - fpr) * (ptpr - tpr) / 2
        #print fpr, tpr, auc, ccc
        pfpr = fpr
        ptpr = tpr

        of.write("%s\t%s\t%s\t%s\t%s\t%f\n" % (p, fpr, tpr, f1, prec, auc))
    of.write("%s\t0.0\t0.0\t0.0\t1.0\t%f\n" % (max_value, auc))

def plot_single(ddir, species, roc_dir, format = "pdf"):
    of = tempfile.NamedTemporaryFile()
    plotfn = "%s/%s/%s" % (ddir, roc_dir, species)
    fn1 = "%s/%s/%s.full" % (ddir, roc_dir, species)
    fn2 = "%s/%s/%s.naive" % (ddir, roc_dir, species)
    fn3 = "%s/%s/%s.blast" % (ddir, roc_dir, species)
    fn4 = "%s/%s/%s.gtg" % (ddir, roc_dir, species)
    fn5 = "%s/%s/%s.blasttree" % (ddir, roc_dir, species)
    fn6 = "%s/%s/%s.gtgtree" % (ddir, roc_dir, species)

    of.write(plot.plot_format(plotfn, format))
    of.write("""
full = read.table("%s")
naive = read.table("%s")
blast = read.table("%s")
gtg = read.table("%s")
btree = read.table("%s")
gtree = read.table("%s")
plot(full[,2], full[,3], type="l", xlab="False positive rate", ylab="True positive rate")
lines(naive[,2], naive[,3], type="l", pch=22, lty=1, col="red")
lines(blast[,2], blast[,3], type="l", pch=22, lty=2, col="blue")
lines(gtg[,2], gtg[,3], type="l", pch=22, lty=2, col="green")
lines(btree[,2], btree[,3], type="l", pch=22, lty=1, col="purple")
lines(gtree[,2], gtree[,3], type="l", pch=22, lty=1, col="orange")
lines(c(0,1), c(0,1))

legend("bottomright", c("Full","Naive","BLAST","GTG", "BlastTree", "GTGTree"), cex=0.8, col=c("black","red","blue","green","purple","orange"), pch=c(21,21,21,21,21,21), lty=c(1,1,2,2,1,1), inset=0.01)
dev.off()
""" % (fn1, fn2, fn3, fn4, fn5, fn6))

    plotfn2 = "%s/%s/%s-F1" % (ddir, roc_dir, species)
    of.write(plot.plot_format(plotfn2, format))
    of.write("""

plot(full[3:nrow(full)-1,1], full[3:nrow(full)-1,4], type="l", ylim=c(0,1), xlab="Normalized score", ylab="")

lines(naive[3:nrow(naive)-1,1], naive[3:nrow(naive)-1,4], col="red", lty=1)
lines(blast[3:nrow(blast)-1,1] / max(blast[,1]), blast[3:nrow(blast)-1,4], col="blue", lty=2)
lines(gtg[3:nrow(gtg)-1,1], gtg[3:nrow(gtg)-1,4], col="green", lty=2)
lines(btree[3:nrow(btree)-1,1], btree[3:nrow(btree)-1,4], col="purple", lty=1)
lines(gtree[3:nrow(gtree)-1,1], gtree[3:nrow(gtree)-1,4], col="orange", lty=1)

legend("bottomleft", c("Full","Naive","BLAST","GTG","BlastTree","GTGTree"), cex=0.8, col=c("black","red","blue","green","purple","orange"), pch=c(21,21,21,21,21,21), lty=c(1,1,2,2,1,1), inset=0.01)
dev.off()
""")

#par(mfrow=c(2,2), oma=c(2, 0, 3, 0), omi=c(0, 0, 0.8, 0))
# plot(naive[3:nrow(naive)-1,1], naive[3:nrow(naive)-1,4], ylim=c(0,1), type="l", xlab="Posterior probability, naive model", ylab="F1 score")
# plot(blast[3:nrow(blast)-1,1] / max(blast[,1]), blast[3:nrow(blast)-1,4], ylim=c(0,1), type="l", xlab="BLAST score", ylab="F1 score")
# plot(gtg[3:nrow(gtg)-1,1], gtg[3:nrow(gtg)-1,4], ylim=c(0,1), type="l", xlab="GTG score", ylab="F1 score")

    plotfn3 = "%s/%s/%s-sens" % (ddir, roc_dir, species)
    of.write(plot.plot_format(plotfn3, format))
    of.write("""
plot(full[1:nrow(full),1], full[1:nrow(full),3], type="l", ylim=c(0,1), xlab="Normalized score", ylab="Sensitivity")

lines(naive[1:nrow(naive),1], naive[1:nrow(naive),3], col="red", lty=1)
lines(blast[1:nrow(blast),1] / max(blast[,1]), blast[1:nrow(blast),3], col="blue", lty=2)
lines(gtg[1:nrow(gtg),1], gtg[1:nrow(gtg),3], col="green", lty=2)
lines(btree[1:nrow(btree),1], btree[1:nrow(btree),3], col="purple", lty=1)
lines(gtree[1:nrow(gtree),1], gtree[1:nrow(gtree),3], col="orange", lty=1)
legend("bottomleft", c("Full","Naive","BLAST","GTG","BlastTree","GTGTree"), cex=0.8, col=c("black","red","blue","green","purple","orange"), pch=c(21,21,21,21,21,21), lty=c(1,1,2,2,1,1), inset=0.01)
dev.off()

""")

    plotfn4 = "%s/%s/%s-spec" % (ddir, roc_dir, species)
    of.write(plot.plot_format(plotfn4, format))
    of.write("""
plot(full[1:nrow(full),1], full[1:nrow(full),5], type="l", ylim=c(0,1), lty=2, xlab="Normalized score", ylab="Specificity")

lines(naive[1:nrow(naive),1], naive[1:nrow(naive),5], col="red", lty=2)
lines(blast[1:nrow(blast),1] / max(blast[,1]), blast[1:nrow(blast),5], col="blue", lty=2)
lines(gtg[1:nrow(gtg),1], gtg[1:nrow(gtg),5], col="green", lty=2)
lines(btree[1:nrow(btree),1], btree[1:nrow(btree),5], col="purple", lty=2)
lines(gtree[1:nrow(gtree),1], gtree[1:nrow(gtree),5], col="orange", lty=2)
legend("bottomleft", c("Full","Naive","BLAST","GTG","BlastTree","GTGTree"), cex=0.8, col=c("black","red","blue","green","purple","orange"), pch=c(21,21,21,21,21,21), lty=c(1,1,1,1,1,1), inset=0.01)
dev.off()

""")

    of.flush()
    subprocess.call("R CMD BATCH %s" % (of.name), shell = True)
    #subprocess.call("R --no-save < %s" % (of.name), shell = True)

def compute_roc_curves(models, scores, all_ecs, roc_dir = common.ROC_CURVE_DIR):
    try:
        os.mkdir("%s/%s" % (ddir, common.ROC_CURVE_DIR))
    except:
        pass
    for species in models:
        print species
        if species not in scores:
            print "Model species %s not in data - skipping" % (species)
            continue
        of = open("%s/%s/%s.full" % (ddir, roc_dir, species), "w")
        roc_curves = compute_roc_curve(models, scores, all_ecs, of, species, 
                                       "pscore", 1)
        of.close()
        of = open("%s/%s/%s.naive" % (ddir, roc_dir, species), "w")
        roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                       "npscore", 1)
        of.close()
        of = open("%s/%s/%s.blast" % (ddir, roc_dir, species), "w")
        roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                       "bscore", None)
        of.close()
        of = open("%s/%s/%s.gtg" % (ddir, roc_dir, species), "w")
        roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                       "gscore", 1)
        of.close()
        of = open("%s/%s/%s.blasttree" % (ddir, roc_dir, species), "w")
        roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                       "btscore", None)
        of.close()
        of = open("%s/%s/%s.gtgtree" % (ddir, roc_dir, species), "w")
        roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                       "gtscore", 1)
        of.close()

        plot_single(ddir, species, roc_dir, "pdf")
        plot_single(ddir, species, roc_dir, "png")

    subprocess.call("rm *.Rout", shell = True)

def roc(ddir, mdir, remove_partial = FILTER_PARTIAL_ECS):
    print "Reading models..."
    models, model_ecs = common.read_models(mdir, remove_partial)
    print "Reading reaction scores..."
    scores, data_ecs = common.read_scores(ddir, True)
    if remove_partial:
        data_ecs = filter_partial_ecs(data_ecs)

    all_ecs = model_ecs.union(data_ecs)
    print "Computing ROC curves..."
    compute_roc_curves(models, scores, all_ecs)

def roc_one(ddir, modelfn, species, outdir, remove_partial = FILTER_PARTIAL_ECS):
    roc_dir = "%s/%s" % (ddir, outdir)
    try:
        os.mkdir(roc_dir)
    except:
        pass

    model = common.read_model(open(modelfn), remove_partial)
    scores, data_ecs = common.read_scores(ddir, True)
    if remove_partial:
        data_ecs = filter_partial_ecs(data_ecs)
    all_ecs = model.union(data_ecs)

    models = {species : model}
    of = open("%s/%s.full" % (roc_dir, species), "w")
    roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                   "pscore", 1, remove_partial)
    of.close()
    of = open("%s/%s.naive" % (roc_dir, species), "w")
    roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                   "npscore", 1, remove_partial)
    of.close()
    of = open("%s/%s.blast" % (roc_dir, species), "w")
    roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                   "bscore", None, remove_partial)
    of.close()
    of = open("%s/%s.gtg" % (roc_dir, species), "w")
    roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                   "gscore", 1, remove_partial)
    of.close()
    of = open("%s/%s.blasttree" % (roc_dir, species), "w")
    roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                   "btscore", 1, remove_partial)
    of.close()
    of = open("%s/%s.gtgtree" % (roc_dir, species), "w")
    roc_curves = compute_roc_curve(models, scores, all_ecs, of, species,
                                   "gtscore", 1, remove_partial)
    of.close()
    plot_single(ddir, species, outdir, "pdf")
    plot_single(ddir, species, outdir, "png")

    subprocess.call("rm *.Rout", shell = True)   

if __name__ == "__main__":
    if len(sys.argv) == 3:
        ddir = sys.argv[1] # project dir
        mdir = sys.argv[2] # model dir
        roc(ddir, mdir)
    elif len(sys.argv) == 6:
        ddir = sys.argv[1] # project dir
        modelfn = sys.argv[2] # model file
        species = sys.argv[3]  # species in reaction score data
        outdir = sys.argv[4]  # output dir
        remove_partial = sys.argv[5] == "yes" # remove partial ECs
        roc_one(ddir, modelfn, species, outdir, remove_partial)
    else:
        print """Compute ROC curves for reaction scores against model files.

Single ROC curve:
   %s project-dir model-file species out-dir remove-partial

Multiple ROC curves:
   %s project-dir model-dir
""" % (sys.argv[0], sys.argv[0])

