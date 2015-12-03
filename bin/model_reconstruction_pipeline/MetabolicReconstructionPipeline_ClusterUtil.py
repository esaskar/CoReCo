#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""


import os
import sys
import traceback
import NGS_Util
import NGS_Blast
import time

sys.path.append("..")
import ScriptsDir


class MetabolicReconstructionPipeline_ClusterUtil:

    orgListFile         = ""  #Name of file containing Organisms List                  #need to be initialized     
    orgFastaDir         = ""  #Directory path containing organisms fasta sequences     #need to be initialized     
    
    orgGTGBlastResDir   = ""  #Directory to store blast outputs                        #need to be initialized     
    GTGFungiBestHitsDir = ""  #Directory to store joint blast rectified outputs        #need to be initialized
    GTGKNNDir           = ""  #Directory to store joint blast rectified outputs        #need to be initialized
    
    orgBlastResDir      = ""  #Directory to store blast outputs                        #need to be initialized     
    jointBlastDir       = ""  #Directory to store joint blast rectified outputs        #need to be initializ

    orgIPRScanDir             = ""
    InterProScan_EC_RAW_results = ""

    
    
    moveToDir = ""
    moveToDir_orgGTGBlastResDir   = ""  #Directory to store blast outputs                        #need to be initialized     
    moveToDir_GTGFungiBestHitsDir = ""  #Directory to store joint blast rectified outputs        #need to be initialized
    moveToDir_GTGKNNDir           = ""  #Directory to store joint blast rectified outputs        #need to be initialized
    
    moveToDir_orgBlastResDir      = ""  #Directory to store blast outputs                        #need to be initialized     
    moveToDir_jointBlastDir       = ""  #Directory to store joint blast rectified outputs        #need to be initializ

    moveToDir_orgIPRScanDir             = ""
    moveToDir_InterProScan_EC_RAW_results = ""

    
    def initialize(self, orgListFile, orgFastaDir,  orgBlastResDir, jointBlastDir, orgGTGBlastResDir, GTGBestHitsDir, GTGKNNDir, orgIPRScanDir, InterProScan_EC_RAW_results, moveToDir): #, new_seq_org_list):
    
        try:
            
            self.orgListFile               = orgListFile
            self.orgFastaDir               = orgFastaDir

            self.orgBlastResDir            = orgBlastResDir
            self.jointBlastDir             = jointBlastDir

	    self.orgGTGBlastResDir         = orgGTGBlastResDir
	    self.GTGBestHitsDir            = GTGBestHitsDir
	    self.GTGKNNDir                 = GTGKNNDir

            self.orgIPRScanDir             = orgIPRScanDir
            self.InterProScan_EC_RAW_results = InterProScan_EC_RAW_results



	    self.moveToDir                             = moveToDir 
	    self.moveToDir_orgBlastResDir              = NGS_Util.createDirectoryPath(self.moveToDir, "blast_results")
	    self.moveToDir_jointBlastDir               = NGS_Util.createDirectoryPath(self.moveToDir, "blast_joint_results") 
	    
	    self.moveToDir_orgGTGBlastResDir           = NGS_Util.createDirectoryPath(self.moveToDir, "GTG_blast_results")    
	    self.moveToDir_GTGBestHitsDir              = NGS_Util.createDirectoryPath(self.moveToDir, "GTG_best_hits")  
	    self.moveToDir_GTGKNNDir                   = NGS_Util.createDirectoryPath(self.moveToDir, "GTG_knn")      
	    
	    self.moveToDir_orgIPRScanDir               = NGS_Util.createDirectoryPath(self.moveToDir, "iprscan_results")  
	    self.moveToDir_InterProScan_EC_RAW_results = NGS_Util.createDirectoryPath(self.moveToDir, "iprscan_ec_raw_results") #InterProScan_EC_RAW_results    


        except Exception:
            
            print traceback.print_exc()
    


        
    def moveFile_createLink(self, srcFile, dstFile):
    
        try:
	    print srcFile + "....." + dstFile

	    if ( not os.path.islink(srcFile) ) and (os.path.isfile(srcFile)):

		NGS_Util.moveFile(srcFile,dstFile)
		os.symlink(dstFile,srcFile)

	    else:

		if ( not os.path.exists(srcFile) ) and ( os.path.exists(dstFile)):

		    os.symlink(dstFile,srcFile)

        except Exception:

            print traceback.print_exc()


        
    def moveIPRScanResults(self):
    
        try:
                
            print "moveIPRScanResults"

            orgListFile_fh = open(self.orgListFile)

            NGS_Util.zipDirectory(self.orgIPRScanDir)
	    NGS_Util.moveDirectoryFiles(self.orgIPRScanDir,self.moveToDir_orgIPRScanDir)

            for line in orgListFile_fh:

                organismNameID, organismName = line.strip().split()

		organism_raw_final        = NGS_Util.createFilePath(self.InterProScan_EC_RAW_results, organismName + ".faa.raw.txt")
		moveto_organism_raw_final = NGS_Util.createFilePath(self.moveToDir_InterProScan_EC_RAW_results, organismName + ".faa.raw.txt")

		self.moveFile_createLink(organism_raw_final,moveto_organism_raw_final)


		organism_IPR_final        = NGS_Util.createFilePath(self.InterProScan_EC_RAW_results, organismName + ".faa.IPR.final.txt")
		moveto_organism_IPR_final = NGS_Util.createFilePath(self.moveToDir_InterProScan_EC_RAW_results, organismName + ".faa.IPR.final.txt")
		
		self.moveFile_createLink(organism_IPR_final,moveto_organism_IPR_final)
		    

            orgListFile_fh.close() 
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""



    def moveGTGResults(self):
    
        try:
                
            print "moveGTGResults"

            orgListFile_fh = open(self.orgListFile)


            NGS_Util.zipDirectory(self.orgGTGBlastResDir)

	    NGS_Util.moveDirectoryFiles(self.orgGTGBlastResDir,self.moveToDir_orgGTGBlastResDir)

            NGS_Util.zipDirectory(self.GTGBestHitsDir)

	    NGS_Util.moveDirectoryFiles(self.GTGBestHitsDir,self.moveToDir_GTGBestHitsDir)
	    
	    
	    for line in orgListFile_fh:

                organismNameID, organismName = line.strip().split()

		org_gtg_knn_final        = NGS_Util.createFilePath(self.GTGKNNDir, organismNameID + ".gtg.knn")
		moveto_org_gtg_knn_final = NGS_Util.createFilePath(self.moveToDir_GTGKNNDir, organismNameID + ".gtg.knn")

		self.moveFile_createLink(org_gtg_knn_final,moveto_org_gtg_knn_final)
		    

            orgListFile_fh.close() 
     
        except Exception:
            
            print traceback.print_exc()
            

    def moveBLASTResults(self):
    
        try:
                
            print "moveBLASTResults"

            orgListFile_fh = open(self.orgListFile)

            NGS_Util.zipDirectory(self.orgBlastResDir)
	    
            NGS_Util.moveDirectoryFiles(self.orgBlastResDir,self.moveToDir_orgBlastResDir)


	    for line in orgListFile_fh:

                organismNameID, organismName = line.strip().split()

		orgRectifyBlast = NGS_Util.createFilePath(self.jointBlastDir, organismName + ".joint.blast")
		moveto_orgRectifyBlast = NGS_Util.createFilePath(self.moveToDir_jointBlastDir, organismName + ".joint.blast")

		self.moveFile_createLink(orgRectifyBlast,moveto_orgRectifyBlast)
		    

            orgListFile_fh.close() 
     
        except Exception:
            
            print traceback.print_exc()


    def moveResults(self):
    
        try:
                
            print "moveResults"
	    
	    self.moveBLASTResults()
	    self.moveGTGResults()
	    self.moveIPRScanResults()
     
        except Exception:
            
            print traceback.print_exc()
            

    def copySequenceFiles(self, srcDataDir):
    
        try:
                
            print("Copy Fasta Files from %s to %s" %(srcDataDir,self.orgFastaDir))

            orgListFile_fh = open(self.orgListFile)

	    for line in orgListFile_fh:

                organismNameID, organismName = line.strip().split()

 		if not os.path.exists( NGS_Util.createFilePath(self.orgFastaDir, organismName + ".faa") ):

 		    orgFasta = NGS_Util.createFilePath(srcDataDir, organismName + ".faa")
		
		    NGS_Util.copyFile(orgFasta, self.orgFastaDir)		    
                    print("Copied fasta file for %s" % (organismName))
                else:
                    print("\tDoing nothing (files already copied) for %s" % (organismName))
                

            orgListFile_fh.close() 
     
        except Exception:
            
            print traceback.print_exc()
	    
	    
