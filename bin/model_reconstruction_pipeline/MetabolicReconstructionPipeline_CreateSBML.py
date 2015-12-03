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

sys.path.append(ScriptsDir.ReconstructionScripts)


class MetabolicReconstructionPipeline_CreateSBML:

    modelTrainingModelDir        = ""
    keggDataDir                  = ""
    taxonomy                     = ""
    intAaccept                   = 0.03
    intReject                    = 2
    networkReconstructionOrgList = ""
    ec2rxnFile                   = ""
   
    def initialize(self, boundsFile, pathwayFile, rxnNamesFile, modelTrainingModelDir, intAaccept, intReject, keggDataDir, ec2rxnFile, taxonomy, networkReconstructionOrgList):
    
        try:
            self.bounds			      = boundsFile	
            self.pathwayFile                  = pathwayFile
            self.rxnNamesFile                 = rxnNamesFile
            self.modelTrainingModelDir        = modelTrainingModelDir
            self.keggDataDir                  = keggDataDir
            self.ec2rxnFile                   = ec2rxnFile
            self.taxonomy                     = taxonomy
            self.intAaccept                   = intAaccept
            self.intReject                    = intReject
            self.networkReconstructionOrgList = networkReconstructionOrgList
    
        except Exception:
            print traceback.print_exc()

       
    def doCreateSBML(self):
    
        try:

            print "doNetworkReconstruction"

	    curDirectory = os.getcwd()
	    os.chdir(ScriptsDir.ReconstructionScripts)
        
            orgListFile_fh = open(self.networkReconstructionOrgList)

            for orgLine in orgListFile_fh:
                
                organismID = orgLine.strip()
		
		print "Network reconstruction for: " + organismID

                intAccept2 = "%.9f" % (self.intAaccept)
                intReject2 = "%.9f" % (self.intReject)
            
		call = "bash " +  ScriptsDir.ReconstructionScripts_reco_dir_postprocess + " " + self.modelTrainingModelDir + " " + intAccept2 + " " + intReject2 + " " + organismID  + " " + ScriptsDir.ReconstructionScripts +  " " +  self.keggDataDir  + " " + self.taxonomy  + " " + self.bounds + " " + self.ec2rxnFile  + " " + "3" +" " + self.rxnNamesFile + " " + self.pathwayFile

		NGS_Util.executeCall(call)

            orgListFile_fh.close() 
     

	    os.chdir(curDirectory)

        except Exception:
            
            print traceback.print_exc()
            
        return ""
    

