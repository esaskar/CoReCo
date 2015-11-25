#!/usr/bin/env python

import sys, os
sys.path.append("../model_training_scripts/")
import common
import reconstruct



def compute_reaction_scores(ddir, r2ec, ec2r, target_species, remove_partial = False):
    scores, all_ecs = common.read_scores(ddir, True)
    try:
        os.mkdir("%s/%s" % (ddir, common.REACTION_SCORE_DIR))
    except:
        pass
    print "Computing reaction scores..."
    for species in scores:
        if target_species != None and species not in target_species:
            continue
        print "   %s" % (species)

        sfn = "%s/%s/%s" % (ddir, common.REACTION_SCORE_DIR, species)
        o = open(sfn, "w")
        rscores = {}
        rec = {}
        for ec in scores[species]:
            if remove_partial and "-" in ec:
                print "reco-dir.py: Discarding reaction score for partial EC number:", ec
                continue
            score = scores[species][ec]
            if ec in ec2r:
                for r in ec2r[ec]:
                    if r not in rec:
                        rec[r] = set()
                    if r not in rscores or score.pscore > rscores[r].pscore:
                        rscores[r] = score
                        rec[r].add(ec)
            else:
                o.write("?\t%s\t%s\n" % (ec, score))

        rs = rscores.keys()
        rs.sort()
        for r in rs:
            o.write("%s\t%s\t%s\n" % (r, ",".join(rec[r]), rscores[r]))
        o.close()
    return scores.keys()

def do_reconstruction(cdir, kdir, ddir, r2ec, ec2r, param_accept, param_reject, target_species = None, odir_param_labeled = False, remove_partial = False):
    species = compute_reaction_scores(ddir, r2ec, ec2r, target_species, remove_partial)

    if target_species == None:
        target_species = species

    adir = "%s/atommaps" % (kdir)
    sourcesfn = "%s/aux/sources-augmented" % (cdir)
    ubiqfn = "%s/aux/empty" % (cdir)

    try:
        os.mkdir("%s/%s" % (ddir, common.RECO_RESULT_DIR))
    except:
        pass
    
    for spec in target_species:
        if spec not in species:
            sys.stderr.write("Unable to reconstruct \"%s\" - not in reaction scores\n" % (spec))
            continue
        print "Reconstructing %s" % (spec)
        scoresfn = "%s/%s/%s" % (ddir, common.REACTION_SCORE_DIR, spec)
        if odir_param_labeled:
            odir = "%s/%s/%s-%.9f-%.9f" % (ddir, common.RECO_RESULT_DIR, spec, param_accept, param_reject)
        else:
            odir = "%s/%s/%s" % (ddir, common.RECO_RESULT_DIR, spec)
        try:
            os.mkdir(odir)
        except:
            pass
        reconstruct.main(cdir, kdir, adir, sourcesfn, ubiqfn, scoresfn, param_accept, param_reject, odir)

def main():
    cdir = sys.argv[1]  # program dir (with ec-list.txt)

    kdir = sys.argv[2]  # kegg dir

    ddir = sys.argv[3]  # data dir
    param_accept = float(sys.argv[4])
    param_reject = float(sys.argv[5])

    species = None
    if len(sys.argv) > 6:
        species = sys.argv[6].split(",")
    
    if len(sys.argv) > 7 and sys.argv[7] == "yes":
        odir_param_labeled = True
    else:
        odir_param_labeled = False

    if len(sys.argv) > 8 and sys.argv[8] == "yes":
        remove_partial = True
    else:
        remove_partial = False

    print "Reading EC list..."
    ec2r, r2ec = common.read_ec_list(open("%s/%s" % (cdir, common.FILE_EC_MAP)))

    do_reconstruction(cdir, kdir, ddir, r2ec, ec2r, param_accept, param_reject, species, odir_param_labeled, remove_partial)

if __name__ == "__main__":
    main()
