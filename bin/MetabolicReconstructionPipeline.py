#!/usr/bin/env python2
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""

import sys, os
sys.path.append("%s/model_reconstruction_pipeline" % (os.path.abspath(os.path.dirname(__file__))))
import NGS_Util
import NGS_Blast
import ScriptsDir
import ProjectDir
from ProjectDir import *
import Blast
import GTG
import IPRScan
import PreProcess
import ModelTraining
import NetworkReconstruction

from config import config

# projectDir        = ProjectDir.projectDir

# uniprot_fasta     = ProjectDir.uniprot_fasta
# uniprot_dust      = ProjectDir.uniprot_dust
# uniprot_blast_db  = ProjectDir.uniprot_blast_db

# ec_files          = ProjectDir.ec_files
# ec2go             = ProjectDir.ec2go
# uniprot_sprot_dat = ProjectDir.uniprot_sprot_dat

# orgListFile       = ProjectDir.orgListFile
# orgFastaDir       = ProjectDir.orgFastaDir

# orgBlastDBDir     = ProjectDir.orgBlastDBDir
# orgBlastDustDir   = ProjectDir.orgBlastDustDir

# orgBlastResDir    = ProjectDir.orgBlastResDir
# jointBlastDir     = ProjectDir.jointBlastDir

# nrdb40_fasta        = ProjectDir.nrdb40_fasta
# nrdb40_blast_db     = ProjectDir.nrdb40_blast_db
# nrdb40_dust         = ProjectDir.nrdb40_dust

# orgGTGDatabaseDir   = ProjectDir.orgGTGDatabaseDir

# orgGTGBlastResDir   = ProjectDir.orgGTGBlastResDir
# GTGBestHitsDir      = ProjectDir.GTGBestHitsDir
# GTGKNNDir           = ProjectDir.GTGKNNDir
 
# CAA1Dir             = ProjectDir.CAA1Dir
# nids_up             = ProjectDir.nids_up
# seq_org_list        = ProjectDir.seq_org_list

numberNearestHits   = int(config["n_nearest_hits"])
blastEValue	        = float(config["blast_e_value"])

# orgIPRScanDir               = ProjectDir.orgIPRScanDir
# InterProScan_EC_RAW_results = ProjectDir.InterProScan_EC_RAW_results

# phylogeneticTreeFile  = ProjectDir.phylogeneticTreeFile
# modelTrainingDir      = ProjectDir.modelTrainingDir

# modelTrainingModelDir = ProjectDir.modelTrainingModelDir
# keggDataDir           = ProjectDir.keggDataDir
# intAaccept            = ProjectDir.intAaccept
# intReject             = ProjectDir.intReject
# taxonomy              = ProjectDir.taxonomy
# networkReconstructionOrgList = ProjectDir.networkReconstructionOrgList

# keggAtomMapsDataDir   = ProjectDir.keggAtomMapsDataDir
# boundFile			  = ProjectDir.boundsPath
# pathwayFile           = ProjectDir.pathwayFile
# rxnNamesFile          = ProjectDir.rxnNamesFile

# ################################################################################################################################################################################

# numberOfFragments = ProjectDir.numberOfFragments

numberOfFragments = int(config["fragments"])

preProcess = PreProcess.PreProcess()
preProcess.initialize(uniprot_fasta, uniprot_dust, uniprot_blast_db, nrdb40_fasta, nrdb40_dust, nrdb40_blast_db, orgListFile, orgFastaDir, seq_org_list, orgBlastDBDir, orgBlastDustDir, ec_files, uniprot_sprot_dat, keggAtomMapsDataDir)
preProcess.preProcess()

blast = Blast.Blast()
blast.initialize(uniprot_fasta, ec_files, uniprot_blast_db, orgListFile, orgFastaDir , orgBlastDBDir, orgBlastDustDir, orgBlastResDir, jointBlastDir, blastEValue, uniprotDBSize)
blast.numberOfFragments = numberOfFragments
blast.getBlastScore() 

gtg = GTG.GTG()
gtg.initialize(nrdb40_blast_db, ec_files, orgListFile, orgFastaDir,  orgBlastDBDir, orgBlastDustDir, orgGTGBlastResDir, GTGBestHitsDir, GTGKNNDir, CAA1Dir, nids_up, seq_org_list, numberNearestHits, orgGTGDatabaseDir, blastEValue)	
gtg.getGTGScore()

iprscan = IPRScan.IPRScan()
iprscan.initialize(orgListFile, orgFastaDir,  orgIPRScanDir, InterProScan_EC_RAW_results, ec2go, seq_org_list, numberOfFragments)
iprscan.getIPRScanScore()

modelTraining = ModelTraining.ModelTraining()
modelTraining.initialize(orgListFile, jointBlastDir, GTGKNNDir, InterProScan_EC_RAW_results, phylogeneticTreeFile, modelTrainingDir)
modelTraining.runModelTrainingAlgo()

networkReconstruction = NetworkReconstruction.NetworkReconstruction()
networkReconstruction.initialize(boundFile, exchangeFile, modelTrainingModelDir, intAaccept, intReject, keggDataDir, ec2rxnFile, taxonomy, networkReconstructionOrgList)
networkReconstruction.doNetworkReconstruction()

createSBML = CreateSBML.CreateSBML()
createSBML.initialize(boundFile, pathwayFile, rxnNamesFile, modelTrainingModelDir, intAaccept, intReject, keggDataDir, ec2rxnFile, taxonomy, networkReconstructionOrgList)
createSBML.doCreateSBML()
        
