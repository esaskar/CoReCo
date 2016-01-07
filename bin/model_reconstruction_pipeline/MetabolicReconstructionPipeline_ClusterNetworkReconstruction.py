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
import MetabolicReconstructionPipeline_NetworkReconstruction

sys.path.append(ScriptsDir.ReconstructionScripts)


class MetabolicReconstructionPipeline_ClusterNetworkReconstruction(MetabolicReconstructionPipeline_NetworkReconstruction.MetabolicReconstructionPipeline_NetworkReconstruction):

       
    def doNetworkReconstruction(self):
    
        try:

            print "doNetworkReconstruction"

	    curDirectory = os.getcwd()
	    os.chdir(ScriptsDir.ReconstructionScripts)
        
            orgListFile_fh = open(self.networkReconstructionOrgList)

            for orgLine in orgListFile_fh:
                
                organismID = orgLine.strip()
		
		print "Network reconstruction for: " + organismID

		call = "qsub -t 1-1 " +  ScriptsDir.ReconstructionScripts_reco_dir_cluster + " " + self.modelTrainingModelDir + " " + str(self.intAaccept) + " " + str(self.intReject) + " " + organismID  + " " + ScriptsDir.ReconstructionScripts +  " " +  self.keggDataDir  + " " + self.taxonomy  + " " + self.bounds + " " + self.ec2rxnFile
		
		NGS_Util.executeCall(call)

            orgListFile_fh.close() 
     

	    os.chdir(curDirectory)

        except Exception:
            
            print traceback.print_exc()
            
        return ""
    

