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
import MetabolicReconstructionPipeline_CreateSBML

sys.path.append(ScriptsDir.ReconstructionScripts)


class MetabolicReconstructionPipeline_ClusterCreateSBML(MetabolicReconstructionPipeline_CreateSBML.MetabolicReconstructionPipeline_CreateSBML):

       
    def doCreateSBML(self):
    
        try:

            print "doCreateSBML"

	    curDirectory = os.getcwd()
	    os.chdir(ScriptsDir.ReconstructionScripts)
        
            orgListFile_fh = open(self.networkReconstructionOrgList)

            for orgLine in orgListFile_fh:
                
                organismID = orgLine.strip()
		
		print "Network reconstruction for: " + organismID

                intAccept2 = "%.9f" % (self.intAaccept)
                intReject2 = "%.9f" % (self.intReject)
            
		call = "qsub -t 1-1 " +  ScriptsDir.ReconstructionScripts_reco_dir_postprocess_cluster + " " + self.modelTrainingModelDir + " " + intAccept2 + " " + intReject2 + " " + organismID  + " " + ScriptsDir.ReconstructionScripts +  " " +  self.keggDataDir  + " " + self.taxonomy  + " " + self.bounds + " " + self.ec2rxnFile   + " " + " " + self.rxnNamesFile + " " + self.pathwayFile 
		
		NGS_Util.executeCall(call)

            orgListFile_fh.close() 
     

	    os.chdir(curDirectory)

        except Exception:
            
            print traceback.print_exc()
            
        return ""
    

