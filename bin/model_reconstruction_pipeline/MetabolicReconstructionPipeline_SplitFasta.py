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
import ProjectDir
import glob
from os import listdir


class MetabolicReconstructionPipeline_SplitFasta:


    splitDataDir         = ""
    organismSplitDataDir = ""
    

    def createSplitDataDir(self):
    
        try:
        
            self.splitDataDir = NGS_Util.createDirectoryPath(ProjectDir.projectResultDir, "split_data")

            NGS_Util.createDirectory(self.splitDataDir)
            
            return self.splitDataDir
            
        except Exception:
            
            print traceback.print_exc()
            
        return ""

    
    def createOrganismSplitDataDir(self, organismName):
    
        try:
        
            self.organismSplitDataDir = NGS_Util.createDirectoryPath(self.splitDataDir, organismName)

            NGS_Util.createDirectory(self.organismSplitDataDir)
            
            return self.organismSplitDataDir
            
        except Exception:
            
            print traceback.print_exc()
            
        return ""
    

    def renameSplitFiles(self, organismName):
    
        try:        
            
            counter = 1
            
            for splitFile in glob.glob(self.organismSplitDataDir + "frag*"):
                
                os.rename(splitFile, organismName + "_" + str(counter))
                
                counter += 1
            
            
        except Exception:
            
            print traceback.print_exc()
            
        return ""
    


    def splitOrganismDataFile(self, organismName, organismDataFile, numberOfFragmemts):
        
        try:
            
            if os.path.exists(organismDataFile) and numberOfFragmemts > 1:
                
                self.createSplitDataDir()
                
                self.organismSplitDataDir = self.createOrganismSplitDataDir(organismName)

                if not os.path.exists( NGS_Util.createDirectoryPath(self.organismSplitDataDir, organismName + "_1" ) ): 
                
                    curDirectory = os.getcwd()

                    os.chdir(self.organismSplitDataDir)
                
                    call = ScriptsDir.fastaSplitN + " " + organismDataFile + " " + str(numberOfFragmemts)

                    NGS_Util.executeCall(call)
                
                    self.renameSplitFiles(organismName)

                    os.chdir(curDirectory)
                
            else:
                
                print  "Organism data file " + organismDataFile + " do not exist or number of fragments must be greater than 1 : " + str(numberOfFragmemts)
                
                
        except Exception:
            
            print traceback.print_exc()
            
        return ""
    
