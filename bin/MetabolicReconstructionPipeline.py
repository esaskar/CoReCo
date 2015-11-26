#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""

import os
import sys
import traceback
sys.path.append("model_reconstruction_pipeline")
import NGS_Util
import NGS_Blast
import ScriptsDir
import ProjectDir
import MetabolicReconstructionPipeline_Blast
import MetabolicReconstructionPipeline_GTG
import MetabolicReconstructionPipeline_IPRScan
import MetabolicReconstructionPipeline_PreProcess
import MetabolicReconstructionPipeline_ModelTraining
import MetabolicReconstructionPipeline_NetworkReconstruction


projectDir        = ProjectDir.projectDir

uniprot_fasta     = ProjectDir.uniprot_fasta
uniprot_dust      = ProjectDir.uniprot_dust
uniprot_blast_db  = ProjectDir.uniprot_blast_db

ec_files          = ProjectDir.ec_files
ec2go             = ProjectDir.ec2go
uniprot_sprot_dat = ProjectDir.uniprot_sprot_dat


orgListFile       = ProjectDir.orgListFile
orgFastaDir       = ProjectDir.orgFastaDir

orgBlastDBDir     = ProjectDir.orgBlastDBDir
orgBlastDustDir   = ProjectDir.orgBlastDustDir

orgBlastResDir    = ProjectDir.orgBlastResDir
jointBlastDir     = ProjectDir.jointBlastDir


nrdb40_fasta        = ProjectDir.nrdb40_fasta
nrdb40_blast_db     = ProjectDir.nrdb40_blast_db
nrdb40_dust         = ProjectDir.nrdb40_dust

orgGTGDatabaseDir   = ProjectDir.orgGTGDatabaseDir

orgGTGBlastResDir   = ProjectDir.orgGTGBlastResDir
GTGBestHitsDir      = ProjectDir.GTGBestHitsDir
GTGKNNDir           = ProjectDir.GTGKNNDir
 
CAA1Dir             = ProjectDir.CAA1Dir
nids_up             = ProjectDir.nids_up
seq_org_list        = ProjectDir.seq_org_list

numberNearestHits = 50
blastEValue	    = ProjectDir.blastEValue


orgIPRScanDir               = ProjectDir.orgIPRScanDir
InterProScan_EC_RAW_results = ProjectDir.InterProScan_EC_RAW_results


phylogeneticTreeFile = ProjectDir.phylogeneticTreeFile
modelTrainingDir     = ProjectDir.modelTrainingDir



modelTrainingModelDir = ProjectDir.modelTrainingModelDir
keggDataDir           = ProjectDir.keggDataDir
intAaccept            = ProjectDir.intAaccept
intReject             = ProjectDir.intReject
taxonomy              = ProjectDir.taxonomy
networkReconstructionOrgList = ProjectDir.networkReconstructionOrgList

keggAtomMapsDataDir = ProjectDir.keggAtomMapsDataDir
boundFile			= ProjectDir.boundsPath
pathwayFile                     = ProjectDir.pathwayFile
rxnNamesFile                    = ProjectDir.rxnNamesFile

################################################################################################################################################################################

numberOfFragments = ProjectDir.numberOfFragments

try:
    
    preProcess = MetabolicReconstructionPipeline_PreProcess.MetabolicReconstructionPipeline_PreProcess()
    preProcess.initialize(uniprot_fasta, uniprot_dust, uniprot_blast_db,  nrdb40_fasta, nrdb40_dust, nrdb40_blast_db, orgListFile, orgFastaDir, seq_org_list, orgBlastDBDir, orgBlastDustDir, ec_files, uniprot_sprot_dat, keggAtomMapsDataDir )
    preProcess.preProcess()


    blast = MetabolicReconstructionPipeline_Blast.MetabolicReconstructionPipeline_Blast()
    blast.initialize(uniprot_fasta_augmented, ec_files, uniprot_blast_db, orgListFile, orgFastaDir , orgBlastDBDir, orgBlastDustDir, orgBlastResDir, jointBlastDir, blastEValue, uniprotDBSize)
    blast.numberOfFragments = numberOfFragments
    blast.getBlastScore() 


    gtg = MetabolicReconstructionPipeline_GTG.MetabolicReconstructionPipeline_GTG()
    gtg.initialize(nrdb40_blast_db, ec_files, orgListFile, orgFastaDir,  orgBlastDBDir, orgBlastDustDir, orgGTGBlastResDir, GTGBestHitsDir, GTGKNNDir, CAA1Dir, nids_up, seq_org_list, numberNearestHits, orgGTGDatabaseDir,blastEValue)	
    gtg.getGTGScore()
    


    iprscan = MetabolicReconstructionPipeline_IPRScan.MetabolicReconstructionPipeline_IPRScan()
    iprscan.initialize(orgListFile, orgFastaDir,  orgIPRScanDir, InterProScan_EC_RAW_results, ec2go, seq_org_list, numberOfFragments)
    iprscan.getIPRScanScore()
    

    modelTraining = MetabolicReconstructionPipeline_ModelTraining.MetabolicReconstructionPipeline_ModelTraining()
    modelTraining.initialize(orgListFile, jointBlastDir, GTGKNNDir, InterProScan_EC_RAW_results, phylogeneticTreeFile, modelTrainingDir)
    modelTraining.runModelTrainingAlgo()

    networkReconstruction = MetabolicReconstructionPipeline_NetworkReconstruction.MetabolicReconstructionPipeline_NetworkReconstruction()
    networkReconstruction.initialize(boundFile, exchangeFile, modelTrainingModelDir, intAaccept, intReject, keggDataDir, ec2rxnFile, taxonomy, networkReconstructionOrgList)
    networkReconstruction.doNetworkReconstruction()
    
    createSBML = MetabolicReconstructionPipeline_CreateSBML.MetabolicReconstructionPipeline_CreateSBML()
    createSBML.initialize(boundFile, pathwayFile, rxnNamesFile, modelTrainingModelDir, intAaccept, intReject, keggDataDir, ec2rxnFile, taxonomy, networkReconstructionOrgList)
    createSBML.doCreateSBML()
        
except Exception:
    
    print traceback.print_exc()












