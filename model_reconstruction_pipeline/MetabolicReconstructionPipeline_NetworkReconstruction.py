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


class MetabolicReconstructionPipeline_NetworkReconstruction:

    modelTrainingModelDir        = ""
    keggDataDir                  = ""
    taxonomy                     = ""
    intAaccept                   = 0.02
    intReject                    = 2
    networkReconstructionOrgList = ""
   
    def initialize(self, modelTrainingModelDir, intAaccept, intReject, keggDataDir, taxonomy, networkReconstructionOrgList):
    
        try:
        
            self.modelTrainingModelDir        = modelTrainingModelDir
            self.keggDataDir                  = keggDataDir
            self.taxonomy                     = taxonomy
            self.intAaccept                   = intAaccept
            self.intReject                    = intReject
            self.networkReconstructionOrgList = networkReconstructionOrgList

    
        except Exception:
            print traceback.print_exc()

       
    def doNetworkReconstruction(self):
    
        try:

            print "doNetworkReconstruction"

	    curDirectory = os.getcwd()
	    os.chdir(ScriptsDir.ReconstructionScripts)
        
            orgListFile_fh = open(self.networkReconstructionOrgList)

            for orgLine in orgListFile_fh:
                
                organismID = orgLine.strip()
		
		print "Network reconstruction for: " + organismID

		call = "bash " +  ScriptsDir.ReconstructionScripts_reco_dir + " " + self.modelTrainingModelDir + " " + str(self.intAaccept) + " " + str(self.intReject) + " " + organismID  + " " + ScriptsDir.ReconstructionScripts +  " " +  self.keggDataDir  + " " + self.taxonomy  
		
		NGS_Util.executeCall(call)

            orgListFile_fh.close() 
     

	    os.chdir(curDirectory)

        except Exception:
            
            print traceback.print_exc()
            
        return ""
    

