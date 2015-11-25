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
import MetabolicReconstructionPipeline_PreProcess_newreactions

import MetabolicReconstructionPipeline_ModelTraining
import MetabolicReconstructionPipeline_NetworkReconstruction

import MetabolicReconstructionPipeline_ClusterBlast
import MetabolicReconstructionPipeline_ClusterIPRScan
import MetabolicReconstructionPipeline_ClusterGTG
import MetabolicReconstructionPipeline_ClusterNetworkReconstruction

import MetabolicReconstructionPipeline_ClusterCreateSBML
import MetabolicReconstructionPipeline_ClusterUtil

projectDir        		= ProjectDir.projectDir
orgGTGDatabaseDir               = ProjectDir.orgGTGDatabaseDir

uniprot_fasta     		= ProjectDir.uniprot_fasta
uniprot_dust      		= ProjectDir.uniprot_dust
uniprot_blast_db  		= ProjectDir.uniprot_blast_db

ec_files          		= ProjectDir.ec_files
ec2go             		= ProjectDir.ec2go
uniprot_sprot_dat 		= ProjectDir.uniprot_sprot_dat
gene2eclist			= ProjectDir.gene2eclist
ec_files_sprot			= ProjectDir.ec_files_sprot

orgListFile      		= ProjectDir.orgListFile
orgFastaDir       		= ProjectDir.orgFastaDir

orgBlastDBDir     		= ProjectDir.orgBlastDBDir
orgBlastDustDir   		= ProjectDir.orgBlastDustDir

orgBlastResDir    		= ProjectDir.orgBlastResDir
jointBlastDir     		= ProjectDir.jointBlastDir


nrdb40_fasta        		= ProjectDir.nrdb40_fasta
nrdb40_blast_db     		= ProjectDir.nrdb40_blast_db
nrdb40_dust         		= ProjectDir.nrdb40_dust

orgGTGDatabaseDir   		= ProjectDir.orgGTGDatabaseDir

orgGTGBlastResDir   		= ProjectDir.orgGTGBlastResDir
GTGBestHitsDir      		= ProjectDir.GTGBestHitsDir
GTGKNNDir           		= ProjectDir.GTGKNNDir
 
CAA1Dir             		= ProjectDir.CAA1Dir
nids_up             		= ProjectDir.nids_up
seq_org_list        		= ProjectDir.seq_org_list

numberNearestHits 		= 50


orgIPRScanDir               	= ProjectDir.orgIPRScanDir
InterProScan_EC_RAW_results 	= ProjectDir.InterProScan_EC_RAW_results


phylogeneticTreeFile 		= ProjectDir.phylogeneticTreeFile
modelTrainingDir     		= ProjectDir.modelTrainingDir



modelTrainingModelDir 		= ProjectDir.modelTrainingModelDir
keggDataDir          		= ProjectDir.keggDataDir
intAaccept            		= ProjectDir.intAaccept
intReject             		= ProjectDir.intReject
taxonomy              		= ProjectDir.taxonomy
networkReconstructionOrgList 	= ProjectDir.networkReconstructionOrgList

keggAtomMapsDataDir 		= ProjectDir.keggAtomMapsDataDir
ec2rxnFile                      = ProjectDir.ec2rxnFile

blastEValue   			= ProjectDir.blastEValue
uniprotDBSize 			= ProjectDir.uniprotDBSize

boundFile			= ProjectDir.boundsPath
exchangeFile		 	= ProjectDir.exchangeReactionsPath
pathwayFile                     = ProjectDir.pathwayFile
rxnNamesFile                    = ProjectDir.rxnNamesFile


################################################################################################################################################################################

numberOfFragments = ProjectDir.numberOfFragments

moveToDir = ProjectDir.projectDirRepository


try:

    call = "export PATH=" + ScriptsDir.ModelTrainingScripts + ":" + projectDir+ ":" + ScriptsDir.Java_PATH +  ":" + ScriptsDir.LibSBML_PATH  +  ":$PATH"
    NGS_Util.executeCall(call)    


    mode = int(sys.argv[1])

    clusterUtil = MetabolicReconstructionPipeline_ClusterUtil.MetabolicReconstructionPipeline_ClusterUtil()
    clusterUtil.initialize(orgListFile, orgFastaDir,  orgBlastResDir, jointBlastDir, orgGTGBlastResDir, GTGBestHitsDir, GTGKNNDir, orgIPRScanDir, InterProScan_EC_RAW_results, moveToDir)

    if mode == 0:

        pathToSeqs = NGS_Util.createFilePath(moveToDir, "data/org_sequence_db/")        
        preProcess = MetabolicReconstructionPipeline_PreProcess.MetabolicReconstructionPipeline_PreProcess()
        preProcess.initialize(uniprot_fasta, uniprot_dust, uniprot_blast_db,  nrdb40_fasta, nrdb40_dust, nrdb40_blast_db, orgListFile, orgFastaDir, seq_org_list, orgBlastDBDir, orgBlastDustDir, ec_files, uniprot_sprot_dat, keggAtomMapsDataDir )
        preProcess.preProcess()

        successStat=True
        print("\nChecking input files that need to be downloaded ...............")
        successStat = successStat and NGS_Util.checkFile(uniprot_fasta)
        successStat = successStat and NGS_Util.checkFile(uniprot_sprot_dat)
        successStat = successStat and NGS_Util.checkFile(ec2go)

        print("\nChecking input files that need to be extracted for GTG ...............")
        successStat = successStat and NGS_Util.checkFile(nrdb40_fasta)
        successStat = successStat and NGS_Util.checkFile(nrdb40_dust)

        print("\n")
        if not successStat:
            raise Exception("Missing required input files")
      
        print("\nChecking Uniprot BLAST input...............")
        successStat = successStat and NGS_Util.checkFile(uniprot_blast_db + ".phd")
        successStat = successStat and NGS_Util.checkFile(uniprot_blast_db + ".psq")
        successStat = successStat and NGS_Util.checkFile(ec_files)
        
        print("\nChecking EC2GO mapping required to process Interpro results...............")
        successStat = successStat and NGS_Util.checkFile(ec2go)
        
        print("\nChecking GTG BLAST input...............")
        successStat = successStat and NGS_Util.checkFile(nrdb40_blast_db)
        successStat = successStat and NGS_Util.checkFile(nrdb40_blast_db)
        
        print("\nChecking Reaction Bag input files")
        successStat = successStat and NGS_Util.checkFile(keggAtomMapsDataDir)
        successStat = successStat and NGS_Util.checkFile(ec2rxnFile)
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(keggDataDir,"kegg-no-general/reaction"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(keggDataDir,"aux/kegg-compounds"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(keggDataDir,"aux/cofactors"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(keggDataDir,"aux/sources"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(keggDataDir,"aux/empty"))
        successStat = successStat and NGS_Util.checkFile(boundFile)
        successStat = successStat and NGS_Util.checkFile(exchangeFile)
        print("\n")
        if not successStat:
            raise Exception("Missing required input files")
                
    if mode > 0 and mode < 3:
        	
        ## GTG

        successStat=True
        print("\nChecking BLAST tools ...............")
        successStat = successStat and NGS_Util.checkFile(ScriptsDir.fastaSplitN)

        print("\nChecking GTG BLAST input...............")
        #successStat = successStat and NGS_Util.checkFile(nrdb40_blast_db + ".phd")
        successStat = successStat and NGS_Util.checkFile(nrdb40_blast_db + ".psq")
        successStat = successStat and NGS_Util.checkFile(ec_files)
        successStat = successStat and NGS_Util.checkFile(orgGTGDatabaseDir)
        successStat = successStat and NGS_Util.checkFile(nids_up)
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(orgGTGDatabaseDir, "nidindex.store"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(orgGTGDatabaseDir, "nidindex.ptr"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(orgGTGDatabaseDir, "vertex.ptr"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(orgGTGDatabaseDir, "CAA1_leafs.store"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(orgGTGDatabaseDir, "CAA1_leafs.ptr"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(orgGTGDatabaseDir, "CAA1_members.store"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(orgGTGDatabaseDir, "CAA1_members.ptr"))
        successStat = successStat and NGS_Util.checkFile(CAA1Dir + ".ptr")
        successStat = successStat and NGS_Util.checkFile(CAA1Dir + ".data")
        print("\nChecking organism sequence input files for GTG BLAST...............")
        successStat = successStat and NGS_Util.checkFile(orgListFile)
        successStat = successStat and NGS_Util.checkFile(orgFastaDir)
        successStat = successStat and NGS_Util.checkFile(seq_org_list)
        print("\n")
        if not successStat:
            raise Exception("Missing required input files")

        gtg = MetabolicReconstructionPipeline_ClusterGTG.MetabolicReconstructionPipeline_ClusterGTG()
        gtg.initialize(nrdb40_blast_db, ec_files, orgListFile, orgFastaDir,  orgBlastDBDir, orgBlastDustDir, orgGTGBlastResDir, GTGBestHitsDir, GTGKNNDir, CAA1Dir, nids_up, seq_org_list, numberNearestHits, orgGTGDatabaseDir,blastEValue)
        gtg.getGTGScore(mode)
    
        ## Interpro

        successStat=True
        print("\nChecking EC2GO required to process inteproscan results...............")
        successStat = successStat and NGS_Util.checkFile(ec2go)
        print("\nChecking organism sequence input files for Interpro scan...............")
        successStat = successStat and NGS_Util.checkFile(orgListFile)
        successStat = successStat and NGS_Util.checkFile(orgFastaDir)
        successStat = successStat and NGS_Util.checkFile(seq_org_list)
        print("\n")
        if not successStat:
            raise Exception("Missing required input files")

        iprscan = MetabolicReconstructionPipeline_ClusterIPRScan.MetabolicReconstructionPipeline_ClusterIPRScan()
        iprscan.initialize(orgListFile, orgFastaDir,  orgIPRScanDir, InterProScan_EC_RAW_results, ec2go, seq_org_list, numberOfFragments)
        iprscan.getIPRScanScore(mode)
 
        ## BLAST against Uniprot

        successStat=True
        print("\nChecking Uniprot BLAST input...............")
        successStat = successStat and NGS_Util.checkFile(uniprot_blast_db + ".phd")
        successStat = successStat and NGS_Util.checkFile(uniprot_blast_db + ".psq")
        successStat = successStat and NGS_Util.checkFile(ec_files)
        print("\nChecking fasta files for 2-way BLAST...............")
        successStat = successStat and NGS_Util.checkFile(uniprot_fasta)
        successStat = successStat and NGS_Util.checkFile(orgListFile)
        successStat = successStat and NGS_Util.checkFile(orgFastaDir)
        print("\n")
        if not successStat:
            raise Exception("Missing required input files")

        blast = MetabolicReconstructionPipeline_ClusterBlast.MetabolicReconstructionPipeline_ClusterBlast()
        blast.initialize(uniprot_fasta, ec_files, uniprot_blast_db, orgListFile, orgFastaDir , orgBlastDBDir, orgBlastDustDir, orgBlastResDir, jointBlastDir, blastEValue, uniprotDBSize)
        blast.numberOfFragments = numberOfFragments
        blast.getBlastScore(mode)       
    

    if mode == 5:
        
        clusterUtil = MetabolicReconstructionPipeline_ClusterUtil.MetabolicReconstructionPipeline_ClusterUtil()
        moveToDirResults=NGS_Util.createFilePath(moveToDir, "results/")
        clusterUtil.initialize(orgListFile, orgFastaDir,  orgBlastResDir, jointBlastDir, orgGTGBlastResDir, GTGBestHitsDir, GTGKNNDir, orgIPRScanDir, InterProScan_EC_RAW_results, moveToDirResults)
        clusterUtil.moveResults()
    
    if mode == 3:
        
        successStat=True
        print("\nChecking for phylogenetic tree required for model building...............")
        successStat = successStat and NGS_Util.checkFile(phylogeneticTreeFile)
        print("\nChecking organism sequence input files for Interpro scan...............")
        successStat = successStat and NGS_Util.checkFile(orgListFile)
        print("\n")
        if not successStat:
            raise Exception("Missing required input files")

        modelTraining = MetabolicReconstructionPipeline_ModelTraining.MetabolicReconstructionPipeline_ModelTraining()
        modelTraining.initialize(orgListFile, jointBlastDir, GTGKNNDir, InterProScan_EC_RAW_results, phylogeneticTreeFile, modelTrainingDir)
        modelTraining.runModelTrainingAlgo()
        
    if mode == 4:
    
        successStat=True
        print("\nChecking Reaction Bag input files")
        successStat = successStat and NGS_Util.checkFile(keggAtomMapsDataDir)
        successStat = successStat and NGS_Util.checkFile(ec2rxnFile)
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(keggDataDir,"kegg-no-general/reaction"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(keggDataDir,"aux/kegg-compounds"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(keggDataDir,"aux/cofactors"))
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(keggDataDir,"aux/empty"))
        successStat = successStat and NGS_Util.checkFile(boundFile)
        successStat = successStat and NGS_Util.checkFile(exchangeFile)
        successStat = successStat and NGS_Util.checkFile(taxonomy)
        print("\n")
        if not successStat:
            raise Exception("Missing required input files")

        networkReconstruction = MetabolicReconstructionPipeline_ClusterNetworkReconstruction.MetabolicReconstructionPipeline_ClusterNetworkReconstruction()
        networkReconstruction.initialize(boundFile, exchangeFile, modelTrainingModelDir, intAaccept, intReject, keggDataDir, ec2rxnFile, taxonomy, networkReconstructionOrgList)
        networkReconstruction.doNetworkReconstruction()
           
    if mode == 6:

        successStat=True
        print("\nChecking Reaction Bag input files")
        successStat = successStat and NGS_Util.checkFile(keggAtomMapsDataDir)
        successStat = successStat and NGS_Util.checkFile(ec2rxnFile)
        successStat = successStat and NGS_Util.checkFile(NGS_Util.createFilePath(keggDataDir,"aux/kegg-compounds"))
        successStat = successStat and NGS_Util.checkFile(boundFile)
        successStat = successStat and NGS_Util.checkFile(exchangeFile)
        successStat = successStat and NGS_Util.checkFile(pathwayFile)
        successStat = successStat and NGS_Util.checkFile(rxnNamesFile)
        successStat = successStat and NGS_Util.checkFile(taxonomy)
        print("\n")
        if not successStat:
            raise Exception("Missing required input files")
    
        createSBML = MetabolicReconstructionPipeline_ClusterCreateSBML.MetabolicReconstructionPipeline_ClusterCreateSBML()  
        createSBML.initialize(boundFile, exchangeFile, pathwayFile, rxnNamesFile, modelTrainingModelDir, intAaccept, intReject, keggDataDir, ec2rxnFile, taxonomy, networkReconstructionOrgList)
        createSBML.doCreateSBML()
           
    
        
except Exception:
    
    print traceback.print_exc()


