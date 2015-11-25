#!/usr/bin/env python

import sys, os

import sys

import common

PARAM_ACCEPT = 1
PARAM_REJECT = 2

PLOT_PARAM = PARAM_ACCEPT
#PLOT_PARAM = PARAM_REJECT

#REMOVE_PARTIAL_ECS = 0

def main(rdir, dir_prefix, ref_model_fn, ofile, remove_partial):
    ref_model_f = open(ref_model_fn)

    of = open(ofile, "w")

    ref_ecs = set()
    for s in ref_model_f:
        ec = s.strip()
        if remove_partial and "-" in ec:
            continue
        ref_ecs.add(ec)

    all_ecs = set()
    f = open(common.FILE_EC_MAP)
    ec2r, r2ec = common.read_ec_list(f)
    all_ecs = set(ec2r.keys())
    if remove_partial:
        all_ecs = map(lambda x: "-" not in x, all_ecs)

    fns = os.listdir(rdir)
    results = {}
    values = []
    for fna in fns:
        if fna.startswith(dir_prefix):
            sys.stdout.write("Processing %s\n" % (fna))
            params = fna[len(dir_prefix):]
            if PLOT_PARAM == PARAM_ACCEPT:
                param = float(params.split("-")[0])  # accept param
            elif PLOT_PARAM == PARAM_REJECT:
                param = float(params.split("-")[1])  # reject param
            else:
                print "Unknown parameter", PLOT_PARAM
                assert(0)

            #print fna, param
            try:
                f = open("%s/%s/%s" % (rdir, fna, common.NETWORK_EC_FILE))
                #f = open("%s/%s/network.ecs.filtered" % (rdir, fna))
            except:
                sys.stderr.write("Unable to open %s\n" % (fna))
                continue
            res_ecs = set()
            for s in f:
                if s.startswith("#"):
                    continue
                ec = s.strip().split()[0]
                if ec == "?":
                    continue
                if remove_partial and "-" in ec:
                    continue
                res_ecs.add(ec)

            tp = len(ref_ecs.intersection(res_ecs)) + 1
            fp = len(res_ecs.difference(ref_ecs)) + 1
            fn = len(ref_ecs.difference(res_ecs)) + 1
            tn = len(all_ecs) - tp - fp - fn + 3

            try:
                tpr = 1.0 * tp / (tp + fn)
                fpn = 1.0 * fp / (fp + tn)

                prec = 1.0 * tp / (tp + fp)
                recall = 1.0 * tp / (tp + fn)

                f1 = 2 * prec * recall / (prec + recall)

                results[param] = (tp, fp, fn, tn, tpr, fpn, f1)

                #values.append((param, tp, fp, fn, tn, tpr, fpn, f1))
                #print fpn, tpr
            except:
                print "Cannot process %s" % (fna)
                raise




    pfpr = ptpr = 1.0
    auc = 0.0

    keys = results.keys()
    keys.sort()
    of.write("#Param TP FP FN TN TPR FPN F1 AUC\n")
    for param in keys:
        tp, fp, fn, tn, tpr, fpr, f1 = results[param]

        auc += (pfpr - fpr) * tpr + (pfpr - fpr) * (ptpr - tpr) / 2
        print fpr, tpr, auc, f1

        of.write("%.10f\t%d\t%d\t%d\t%d\t%f\t%f\t%s\t%f\n" % (param, tp, fp, fn, tn, tpr, fpr, f1, auc))
        pfpr = fpr
        ptpr = tpr

    tpr = 0.0
    fpr = 0.0
    auc += (pfpr - fpr) * tpr + (pfpr - fpr) * (ptpr - tpr) / 2
    f1 = 0.0
    of.write("NA\tNA\tNA\tNA\tNA\t%f\t%f\t%f\t%f\n" % (tpr, fpr, f1, auc))
    #print fpr, tpr, auc, f1

if __name__ == "__main__":
    rdir = sys.argv[1]  # dir with result dirs
    dir_prefix = sys.argv[2] # result dir prefix, e.g., "t"
    ref_model_fn = sys.argv[3] # reference model ECs
    ofile = sys.argv[4]  # output file
    remove_partial = sys.argv[5] == "yes" # remove partial ECs, e.g., "1.2.-.-"
    main(rdir, dir_prefix, ref_model_fn, ofile, remove_partial)


