#!/usr/bin/env python

import sys, os
sys.path.append("../model_training_scripts/")
import common
import reconstruct



def compute_reaction_scores(ddir, r2ec, ec2r, param_accept, param_reject, target_species, odir_param_labeled, remove_partial):
    scores, all_ecs = common.read_scores(ddir, True)
    print "Computing reaction scores..."
    for species in scores:
        if target_species != None and species not in target_species:
            continue
        print "   %s" % (species)

        sfndir = "%s/%s" % (ddir, common.REACTION_SCORE_DIR)


        try:
            os.mkdir(sfndir)
        except:
            pass
        sfn = "%s/%s" % (sfndir, species)
        
        alreadycomputed=False
        # check if reaction scores already extracted
        if os.path.isfile(sfn):
            if os.path.getsize(sfn) > 0:
                alreadycomputed=True

        if not alreadycomputed:
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
                if not alreadycomputed:
                    o.write("?\t%s\t%s\n" % (ec, score))

        rs = rscores.keys()
        rs.sort()
        if not alreadycomputed:
            for r in rs:
                o.write("%s\t%s\t%s\n" % (r, ",".join(rec[r]), rscores[r]))
            o.close()
    return scores.keys()

def do_reconstruction(bfile, cdir, kdir, ddir, r2ec, ec2r, param_accept, param_reject, cnamefile, ecfile, target_species = None, odir_param_labeled = False, remove_partial = False):
    species = compute_reaction_scores(ddir, r2ec, ec2r, param_accept, param_reject, target_species, odir_param_labeled, remove_partial)

    if target_species == None:
        target_species = species

    adir = "%s/atommaps" % (kdir)
    sourcesfn = "%s/aux/sources" % (cdir)
    ubiqfn = "%s/aux/empty" % (cdir)
    
    if odir_param_labeled:
        odir = "%s/%s-%.9f-%.9f" % (ddir, common.RECO_RESULT_DIR, param_accept, param_reject)
    else:
        odir = "%s/%s" % (ddir, common.RECO_RESULT_DIR)
    try:
        os.mkdir(odir)
    except:
        pass

    for spec in target_species:
        if spec not in species:
            sys.stderr.write("Unable to reconstruct \"%s\" - not in reaction scores\n" % (spec))
            continue
        print "Reconstructing %s" % (spec)
        # reactions scores the same regardless of the reconstruction parameter
        scoresfn = "%s/%s/%s" % (ddir, common.REACTION_SCORE_DIR, spec)
        if odir_param_labeled:
            odir = "%s/%s-%.9f-%.9f/%s" % (ddir, common.RECO_RESULT_DIR, param_accept, param_reject, spec)
        else:
            odir = "%s/%s/%s" % (ddir, common.RECO_RESULT_DIR, spec)
        print "Writing results to %s" % odir
        try:
            os.mkdir(odir)
        except:
            pass
        reconstruct.main(bfile, cdir, kdir, adir, sourcesfn, ubiqfn, scoresfn, param_accept, param_reject, odir, cnamefile, ecfile)

def main():   
	
    bfile = sys.argv[1]	# bounds file

    cdir = sys.argv[2]  # program dir (with ec-list.txt)   

    kdir = sys.argv[3]  # kegg dir

    ddir = sys.argv[4]  # data dir
    param_accept = float(sys.argv[5])
    param_reject = float(sys.argv[6])

    ecfile = sys.argv[7] # ec file
    print ecfile

    cnamefile = ("%s/../aux/kegg-compounds" % kdir )

    species = None
    if len(sys.argv) > 8:
        species = sys.argv[8].split(",")
    
    if len(sys.argv) > 9 and sys.argv[9] == "yes":
        odir_param_labeled = True
    else:
        odir_param_labeled = False

    if len(sys.argv) > 10 and sys.argv[10] == "yes":
        remove_partial = True
    else:
        remove_partial = False

    print "Reading EC list..."
    ec2r, r2ec = common.read_ec_list(open(ecfile))

    

    do_reconstruction(bfile, cdir, kdir, ddir, r2ec, ec2r, param_accept, param_reject, cnamefile, ecfile, species, odir_param_labeled, remove_partial)

if __name__ == "__main__":
    main()
