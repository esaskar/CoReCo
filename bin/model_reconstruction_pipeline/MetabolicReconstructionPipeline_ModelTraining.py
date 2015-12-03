#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""


import os
import sys
import traceback
import NGS_Util
import NGS_Blast
sys.path.append("..")
import ScriptsDir

sys.path.append(ScriptsDir.ModelTrainingScripts)

class MetabolicReconstructionPipeline_ModelTraining:


    seq_org_list  = ""  #need to be absolute
    
    jointBlastDir = ""  #need to be absolute
    
    GTGFungiKNNDir = ""  #need to be absolute
    
    treeCPDS       = ""  #TREE_CPD_FILE = "tree.cpds"
    
    phylogeneticTreeFile = ""  #need to be absolute 
    
    modelTrainingDir = ""

    modelTraining_IPR_EC_Dir = ""

    modelTrainingBlastPVDir = ""
    
    modelTraining_EC_Scores_Dir = ""
    
    modelTrainingProbabilityDensityScoreDir = ""
    
    modelTrainingTreeDir = ""
    modelTrainingModelDir  = ""


    def initialize(self, seq_org_list, jointBlastDir, GTGFungiKNNDir, fungi_InterProScan_result, phylogeneticTreeFile, modelTrainingDir):
    
        try:
        
	    self.seq_org_list                            = seq_org_list

	    self.jointBlastDir                           = jointBlastDir
	    self.GTGFungiKNNDir                          = GTGFungiKNNDir
	    self.fungi_InterProScan_result               = fungi_InterProScan_result	    

	    self.phylogeneticTreeFile                    = phylogeneticTreeFile
	    
	    self.modelTrainingDir                        = modelTrainingDir
	    self.modelTraining_IPR_EC_Dir                = NGS_Util.createDirectoryPath(self.modelTrainingDir, "IPR_EC")
	    self.modelTrainingBlastPVDir                 = NGS_Util.createDirectoryPath(self.modelTrainingDir, "BlastPValues")
	    self.modelTraining_EC_Scores_Dir             = NGS_Util.createDirectoryPath(self.modelTrainingDir, "ECScores")    
	    self.modelTrainingProbabilityDensityScoreDir = NGS_Util.createDirectoryPath(self.modelTrainingDir, "ProbabilityDensityScore")	    
	    self.modelTrainingTreeDir                    = NGS_Util.createDirectoryPath(self.modelTrainingDir, "Tree")	    
	    self.modelTrainingModelDir                   = NGS_Util.createDirectoryPath(self.modelTrainingDir, "Model")


	    NGS_Util.createDirectory(self.modelTrainingDir)
	    NGS_Util.createDirectory(self.modelTraining_IPR_EC_Dir)
	    NGS_Util.createDirectory(self.modelTrainingBlastPVDir)
	    NGS_Util.createDirectory(self.modelTraining_EC_Scores_Dir)
	    NGS_Util.createDirectory(self.modelTrainingProbabilityDensityScoreDir)
	    NGS_Util.createDirectory(self.modelTrainingTreeDir)
	    NGS_Util.createDirectory(self.modelTrainingModelDir)
	    
	    
	    if (os.path.exists(self.phylogeneticTreeFile)):
		
		NGS_Util.copyFile( self.phylogeneticTreeFile,NGS_Util.createFilePath(self.modelTrainingTreeDir,"tree")   )
		self.phylogeneticTreeFile = NGS_Util.createFilePath(self.modelTrainingTreeDir,"tree")


	    self.treeCPDS = NGS_Util.createFilePath(self.modelTrainingTreeDir,"tree.cpds")
    
        except Exception:
            print traceback.print_exc()
    
    
    def extract_IPR_EC_Values(self):
       
        try:
        
            print "extract_IPR_EC_Values" 
	    
	    call = "python " + ScriptsDir.ModelTrainingScripts_extract_ecs_from_iprscan + " " + self.fungi_InterProScan_result  + " " + self.modelTraining_IPR_EC_Dir
	    
	    NGS_Util.executeCall(call)
	    

        except Exception:
            
            print traceback.print_exc()
        

    def computeBlastPvalues(self):
       
        try:
        
            print "computeBlastPvalues"

	    call = "python " + ScriptsDir.ModelTrainingScripts_computeBlastPvalues + " " + self.jointBlastDir + " " + self.modelTrainingBlastPVDir
	    
	    NGS_Util.executeCall(call)
	    

        except Exception:
            
            print traceback.print_exc()
        


    def computeMergedScores(self):
       
        try:
        
            print "computeMergedScores"

	    call = "python " + ScriptsDir.ModelTrainingScripts_computeMergedScores + " " + self.seq_org_list + " " + self.modelTrainingBlastPVDir + " " + self.GTGFungiKNNDir  + " " + self.modelTraining_EC_Scores_Dir 

	    
	    NGS_Util.executeCall(call)
	    
        except Exception:
            
            print traceback.print_exc()
       

    def computeECscores(self):
       
        try:
        
            print "computeECscores"

	    call = "python " + ScriptsDir.ModelTrainingScripts_computeECscores + " " + self.modelTraining_EC_Scores_Dir 
	    
	    NGS_Util.executeCall(call)
	    
        except Exception:
            
            print traceback.print_exc()



    def computeProbabilityDensityScores(self):
       
        try:
        
            print "computeProbabilityDensityScores"

	    call = "python " + ScriptsDir.ModelTrainingScripts_estimate_cpds_wrapper + " " + self.seq_org_list + " " + self.modelTraining_EC_Scores_Dir  + " " + self.modelTraining_IPR_EC_Dir  + " " + self.modelTrainingProbabilityDensityScoreDir
	    
	    NGS_Util.executeCall(call)
	    
        except Exception:
            
            print traceback.print_exc()



    def combineProbabilityDensityScores(self):
       
        try:
        
            print "combineProbabilityDensityScores"

	    call = "python " + ScriptsDir.ModelTrainingScripts_combined_density + " " + self.modelTrainingProbabilityDensityScoreDir
	    
	    NGS_Util.executeCall(call)
	    
        except Exception:
            
            print traceback.print_exc()


	    
    def computeTreeProbabilityDensityScore(self):
       
        try:
        
            print "computeTreeProbabilityDensityScore"

	    call = "python " + ScriptsDir.ModelTrainingScripts_estimate_mutation_probability + " " + self.modelTraining_IPR_EC_Dir  + " " + self.phylogeneticTreeFile   + " " +  self.modelTrainingTreeDir
	    NGS_Util.executeCall(call)

	    self.treeCPDS = NGS_Util.createFilePath(self.modelTrainingTreeDir,"tree.cpds")

        except Exception:
            
            print traceback.print_exc()


    def runBashScript(self):
       
        try:
        
            print "runBashScript"


	    call = "bash " + ScriptsDir.ModelTrainingScripts_run_job + " " + self.modelTrainingDir  + " " + self.modelTrainingModelDir  + " " + self.modelTraining_EC_Scores_Dir   + " " + self.modelTrainingProbabilityDensityScoreDir+"all" + " " + self.seq_org_list +  " " + self.phylogeneticTreeFile + " " + self.modelTraining_IPR_EC_Dir  + " " + self.treeCPDS
	    
	    NGS_Util.executeCall(call)

    
        except Exception:
            
            print traceback.print_exc()



    def runModelTrainingAlgo(self):
       
        try:
        
            print "runModelTrainingAlgo"

	    curDirectory = os.getcwd()

	    os.chdir(ScriptsDir.ModelTrainingScripts)

	    self.extract_IPR_EC_Values()
	    
	    self.computeBlastPvalues()    
	    
	    self.computeMergedScores()

	    self.computeECscores()
	    
	    self.computeProbabilityDensityScores()	

	    self.combineProbabilityDensityScores()

	    self.computeTreeProbabilityDensityScore()

	    self.runBashScript()

	    NGS_Util.executeCall("rm *.Rout")
	    NGS_Util.executeCall("rm *.RData")
	    
	    os.chdir(curDirectory)

        except Exception:
            
            print traceback.print_exc()



