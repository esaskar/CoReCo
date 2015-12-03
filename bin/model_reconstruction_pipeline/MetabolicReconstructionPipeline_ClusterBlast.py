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
import MetabolicReconstructionPipeline_SplitFasta
import MetabolicReconstructionPipeline_Blast

sys.path.append("..")
import ScriptsDir


class MetabolicReconstructionPipeline_ClusterBlast(MetabolicReconstructionPipeline_Blast.MetabolicReconstructionPipeline_Blast):
    
    numberOfFragments = 10
    
    def run_Org_vs_Uniprot_ClusterBlast(self, organismName):
    
        try:
              
            splitFasta = MetabolicReconstructionPipeline_SplitFasta.MetabolicReconstructionPipeline_SplitFasta()
           
            ###############################################################################################################################################

            org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName+".faa")

            splitFasta.splitOrganismDataFile(organismName, org_fasta, self.numberOfFragments)

            ###################################################################################################################################
            
            clusterArrayCall =  "qsub -t 1-" + str(self.numberOfFragments) + ":1 " + ScriptsDir.ClusterBlast
            
            blastP = NGS_Util.createFilePath(ScriptsDir.BlastDir,"blastp")
            
            outfmt = str(6)
            
            splitFile = splitFasta.organismSplitDataDir + organismName 

            org_vs_UniprotBlastDB =  NGS_Util.createFilePath(self.orgBlastResDir, organismName + "-vs-up" )
                

            call = clusterArrayCall + " "  + blastP  + " "  + self.uniprot_blast_db + " "  + splitFile + " " + outfmt + " " + org_vs_UniprotBlastDB  + " " + str(self.blastEValue)
                

            NGS_Util.executeCall(call)



        except Exception:
            
            print traceback.print_exc()
            




    def concatenate_Org_vs_Uniprot_ClusterBlast_results(self, organismName):
    
        try:
            
            clusterProcessing = True
            
            for fragment in range(self.numberOfFragments):
                                           
                org_vs_UniprotBlastDB =  NGS_Util.createFilePath(self.orgBlastResDir, organismName + "-vs-up_" + str(fragment+1) + ".blast" )
                
                if not os.path.exists(org_vs_UniprotBlastDB):
                    clusterProcessing = False
                    break
                    
                   
            if (clusterProcessing):
                                   
                org_vs_UniprotBlastDB =  NGS_Util.createFilePath(self.orgBlastResDir, organismName + "-vs-up.blast" )
    
                call = "cat " + NGS_Util.createFilePath(self.orgBlastResDir, organismName + "-vs-up_*") + " > " + org_vs_UniprotBlastDB
    
                NGS_Util.executeCall(call)
                
                
                return org_vs_UniprotBlastDB
            
            else:
                print organismName + "-vs-Uniprot BLAST incomplete"
            
    
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def run_Uniprot_vs_Org_ClusterBlast(self, organismName):
    
        try:

                
            splitFasta = MetabolicReconstructionPipeline_SplitFasta.MetabolicReconstructionPipeline_SplitFasta()
           
            ###################################################################################################################################

            splitFasta.splitOrganismDataFile("uniprot", self.uniprot_fasta, self.numberOfFragments)
            
            org_blast_db = NGS_Util.createFilePath(self.orgBlastDBDir, organismName)

            ###################################################################################################################################

            clusterArrayCall =  "qsub -t 1-" + str(self.numberOfFragments) + ":1 " + ScriptsDir.ClusterBlast

           
            blastP = NGS_Util.createFilePath(ScriptsDir.BlastDir,"blastp")
            
            outfmt = str(6)
                        
            
            splitFile = splitFasta.organismSplitDataDir + "uniprot"

            Uniprot_vs_orgBlastDB = NGS_Util.createFilePath(self.orgBlastResDir, "up-vs-" + organismName )

            
            call = clusterArrayCall + " "  + blastP  + " "  + org_blast_db + " "  + splitFile + " " + outfmt + " " + Uniprot_vs_orgBlastDB  + " " + str(self.blastEValue) # +  " " + str(self.uniprotDBSize)


            NGS_Util.executeCall(call)
            
    
        except Exception:
            
            print traceback.print_exc()
            



    def concatenate_Uniprot_vs_Org_ClusterBlast_results(self, organismName):
    
        try:

          
            clusterProcessing = True
            
            for fragment in range(self.numberOfFragments):
                                           
                Uniprot_vs_orgBlastDB = NGS_Util.createFilePath(self.orgBlastResDir, "up-vs-" + organismName + "_" + str(fragment+1) + ".blast")
                
                if not os.path.exists(Uniprot_vs_orgBlastDB):
                    clusterProcessing = False
                    break
                    
                   
            if (clusterProcessing):
                                                           
                Uniprot_vs_orgBlastDB = NGS_Util.createFilePath(self.orgBlastResDir, "up-vs-" + organismName + ".blast")
    
                call = "cat " + NGS_Util.createFilePath(self.orgBlastResDir, "up-vs-" + organismName + "*") + " > " + Uniprot_vs_orgBlastDB
                
                NGS_Util.executeCall(call)
                
                
                return Uniprot_vs_orgBlastDB

            else:
                print "Uniprot-vs-" + organismName + " blast incomplete for:" +  organismName
                          
 
        except Exception:
            
            print traceback.print_exc()
            
        return ""



    def getBlastScore(self, mode):
    
            try:

                print "getBlastScores"

                orgListFile_fh = open(self.orgListFile)
    
                for line in orgListFile_fh:
                    
                    organismNameID, organismName = line.strip().split()
                    
                    orgJointBlast   = NGS_Util.createFilePath(self.orgBlastResDir, organismName + ".joint.blast")
                    orgRectifyBlast = NGS_Util.createFilePath(self.jointBlastDir, organismName + ".joint.blast")
                    
                    print "getBlastScore:" + organismName
                    if not os.path.exists(orgRectifyBlast):
                        
                        if os.path.exists(orgJointBlast):

                            orgRectifyBlast = self.rectifyBlast(organismName, orgJointBlast)
                        
                        else:

                            if (mode == 1):                                                        
                                
                                org_blast_db = self.makeBlastDB(organismName)
                                                       
                                self.run_Org_vs_Uniprot_ClusterBlast(organismName)
                                
                                time.sleep(1800) #wait for 15 minutes
                                    
                                self.run_Uniprot_vs_Org_ClusterBlast(organismName)

                                time.sleep(2400) #wait for 20 minutes

                            elif (mode == 2):
                                
                                org_vs_UniprotBlastDB = self.concatenate_Org_vs_Uniprot_ClusterBlast_results(organismName)
                                Uniprot_vs_orgBlastDB = self.concatenate_Uniprot_vs_Org_ClusterBlast_results(organismName)
                            
                                if (org_vs_UniprotBlastDB != "" and Uniprot_vs_orgBlastDB != ""):
                                    orgJointBlast = self.combineBlast(organismName, org_vs_UniprotBlastDB, Uniprot_vs_orgBlastDB)
                                
                                    if (orgJointBlast != ""):
                                        orgRectifyBlast = self.rectifyBlast(organismName, orgJointBlast)
                              
                orgListFile_fh.close() 
         
            except Exception:
                
                print traceback.print_exc()
                
            return ""
