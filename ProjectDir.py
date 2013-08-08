#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""

import sys
sys.path.append("model_reconstruction_pipeline")
import NGS_Util
import ScriptsDir

projectDir       = ScriptsDir.projectDir  #Needs to be set the user#


projectBinDir    = NGS_Util.createDirectoryPath(projectDir,"bin")
projectDataDir   = NGS_Util.createDirectoryPath(projectDir,"data")
projectResultDir = NGS_Util.createDirectoryPath(projectDir,"results")



##########################################################################################     Blast Toolkit      ##########################################################################################

orgBlastDBDir     = ScriptsDir.BlastDBDir
orgBlastDustDir   = ScriptsDir.BlastDustDir


##########################################################################################     Uniprot Data Files ##########################################################################################

uniprot_fasta     = NGS_Util.createFilePath(projectDataDir, "Uniprot_EC_GO_Data/uniprot_sprot.fasta")
uniprot_sprot_dat = NGS_Util.createFilePath(projectDataDir, "Uniprot_EC_GO_Data/uniprot_sprot.dat")     

uniprot_dust     = NGS_Util.createFilePath(orgBlastDustDir, "uniprot_dust.asnb")
uniprot_blast_db = NGS_Util.createFilePath(orgBlastDBDir, "uniprot")
##########################################################################################     EC Go Data Files   ##########################################################################################

ec_files          = NGS_Util.createFilePath(projectDataDir, "Uniprot_EC_GO_Data/ec_files.txt")
ec2go             = NGS_Util.createFilePath(projectDataDir, "Uniprot_EC_GO_Data/ec2go.txt")      

##########################################################################################     Sequence Data      ##########################################################################################

orgListFile       = NGS_Util.createFilePath(projectDataDir, "org_list.backup")
orgFastaDir       = NGS_Util.createDirectoryPath(projectDataDir, "org_sequence_db")      #"50fungi+newPichia_fasta_seq")      
seq_org_list      = NGS_Util.createFilePath(projectDataDir, "seq_org_list.txt")

taxonomy          = NGS_Util.createFilePath(projectDataDir, "taxonomy")

nrdb40_fasta      = NGS_Util.createFilePath(projectDataDir, "GTG_database/nrdb40_v2.fasta")
orgGTGDatabaseDir = NGS_Util.createDirectoryPath(projectDataDir, "GTG_database")
CAA1Dir           = NGS_Util.createFilePath(projectDataDir, "GTG_database/index/CAA1.index")
nids_up           = NGS_Util.createFilePath(projectDataDir, "GTG_database/nids.up")     



nrdb40_dust     = NGS_Util.createFilePath(orgBlastDustDir, "nrdb40_dust.asnb")
nrdb40_blast_db = NGS_Util.createFilePath(orgBlastDBDir, "nrdb40")

##########################################################################################     Sequence Blast Output  ######################################################################################

orgBlastResDir    = NGS_Util.createDirectoryPath(projectResultDir, "blast_results")
jointBlastDir     = NGS_Util.createDirectoryPath(projectResultDir, "blast_joint_results") 

NGS_Util.createDirectory(orgBlastResDir)
NGS_Util.createDirectory(jointBlastDir)

##########################################################################################     Sequence GTG Output  ######################################################################################

orgGTGBlastResDir   = NGS_Util.createDirectoryPath(projectResultDir, "GTG_blast_results")    
GTGBestHitsDir = NGS_Util.createDirectoryPath(projectResultDir, "GTG_best_hits")  
GTGKNNDir      = NGS_Util.createDirectoryPath(projectResultDir, "GTG_knn")      


NGS_Util.createDirectory(orgGTGBlastResDir)
NGS_Util.createDirectory(GTGBestHitsDir)
NGS_Util.createDirectory(GTGKNNDir)


##########################################################################################     Sequence IPR Output  ######################################################################################


orgIPRScanDir               = NGS_Util.createDirectoryPath(projectResultDir, "iprscan_results")  
InterProScan_EC_RAW_results = NGS_Util.createDirectoryPath(projectResultDir, "iprscan_ec_raw_results") #fungi_InterProScan_result

NGS_Util.createDirectory(orgIPRScanDir)
NGS_Util.createDirectory(InterProScan_EC_RAW_results)

##########################################################################################     Model Training Output  ######################################################################################

#phylogeneticTreeFile = "/mnt/nfs/coreco_final_1/data/tree.txt"
phylogeneticTreeFile = NGS_Util.createFilePath(projectDataDir, "tree.txt")


modelTrainingDir          = NGS_Util.createDirectoryPath(projectResultDir,"ModelTraining")
NGS_Util.createDirectory(modelTrainingDir)

##########################################################################################     Model Rencontruction  ######################################################################################


modelTrainingModelDir        = NGS_Util.createFilePath(modelTrainingDir,"Model")
intAaccept                   = 0.02
intReject                    = 2
keggDataDir                  = NGS_Util.createFilePath(projectDataDir, "Kegg")  
networkReconstructionOrgList = NGS_Util.createFilePath(projectDataDir, "network_reconstruction_org_list.txt") # 4 letter org name file

keggAtomMapsDataDir          = NGS_Util.createDirectoryPath(projectDataDir, "Kegg/kegg-no-general/atommaps") 

