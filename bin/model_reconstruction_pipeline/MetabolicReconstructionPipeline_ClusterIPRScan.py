#!/usr/bin/env python

#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""


import os
import sys
import traceback
import time
import glob
import NGS_Util
import NGS_Blast
import MetabolicReconstructionPipeline_SplitFasta
import MetabolicReconstructionPipeline_IPRScan

sys.path.append("..")
import ScriptsDir


class MetabolicReconstructionPipeline_ClusterIPRScan(MetabolicReconstructionPipeline_IPRScan.MetabolicReconstructionPipeline_IPRScan):


    def runClusterIPRScan(self, organismName):
    
        try:
           
            splitFasta = MetabolicReconstructionPipeline_SplitFasta.MetabolicReconstructionPipeline_SplitFasta()
           
            #numberOfFragments = 10
                       
            ###############################################################################################################################################

            org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName+".faa")

            splitFasta.splitOrganismDataFile(organismName, org_fasta, self.numberOfFragments)

            ###################################################################################################################################
            
            clusterArrayCall =  "qsub -t 1-" + str(self.numberOfFragments) + ":1 " +  ScriptsDir.ClusterIprscan
            
            iprscan = NGS_Util.createFilePath(ScriptsDir.IprscanDir,"interproscan.sh ")
            
            splitFile = splitFasta.organismSplitDataDir + organismName
            
            ipr_raw_file_split = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + "_split" )
                           
            call = clusterArrayCall + " "  + iprscan  + " "  + splitFile + " " + ipr_raw_file_split
            
            NGS_Util.executeCall(call)
            
    
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def mergeXML(self,organismName):
    
        try:

            isFirst = True
            
            ipr_xml_file = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + ".xml")
            ipr_xml_file_fh = open(ipr_xml_file,"w")
            
            for srcFile in glob.glob( NGS_Util.createFilePath(self.orgIPRScanDir, organismName + "_split_*") ):
        
                ipr_XML_split_fh = open(srcFile)
    
                for line in ipr_XML_split_fh:
                                        
                    if line.startswith("<?xml version=") or line.startswith("<protein-matches"):

                        if isFirst:
                            
                            ipr_xml_file_fh.write(line)
                            
                    elif not line.startswith("</protein-matches>"):
                        isFirst = False                         
                        ipr_xml_file_fh.write(line)


                ipr_XML_split_fh.close()


            ipr_xml_file_fh.write("</protein-matches>")
            ipr_xml_file_fh.close()
            
            
            return ipr_xml_file

     
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def concatenate_ClusterIPRScan_results(self, organismName):
    
        try:          
          
            #numberOfFragments = 10
            
            clusterProcessing = True

            for fragment in range(self.numberOfFragments):
                                           
                ipr_raw_file_split = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + "_split_" +  str(fragment+1) + ".xml")
                
                if not os.path.exists(ipr_raw_file_split):                   
                    clusterProcessing = False
                    break
            
            if clusterProcessing:
                
                ipr_xml_file = self.mergeXML(organismName)
                return ipr_xml_file
                
            
            else:
                print "Interpro incomplete for: " +  organismName
                
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def xmlIPRScanToRAWOutput(self, organismName, ipr_xml_file):
    
        try:
        
            print "xmlIPRScanToRAWOutput: " + organismName
            
            ipr_raw_file = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + ".ipr.raw")
            
            self.ngsIPRScan.convert_iprscan5_xml_raw(ipr_xml_file, ipr_raw_file)
            
            return ipr_raw_file

        except Exception:
            
            print traceback.print_exc()
        
        return ""

    def getIPRScanScore(self, mode):
    
        try:
        
                
            print "getIPRScanScore"

            orgListFile_fh = open(self.orgListFile)

            for line in orgListFile_fh:
                
                organismNameID, organismName = line.strip().split()
                
                organism_IPR_final = NGS_Util.createFilePath(self.fungi_InterProScan_result, organismName + ".faa.IPR.final.txt")

               
                if not os.path.exists(organism_IPR_final):
                    
                    print "getIPRScanScore : " + organismName
                    
                    if mode == 1:

                        self.runClusterIPRScan(organismName)
                        
                        time.sleep(21600) # sleep for 6 hrs
                        
                    elif mode == 2:
                        
                        ipr_xml_file = self.concatenate_ClusterIPRScan_results(organismName)                        
                        ipr_raw_file = self.xmlIPRScanToRAWOutput(organismName, ipr_xml_file)
    
                        organism_ipr2go = self.extract_ipr2go_based_on_xml(organismName, ipr_xml_file)
                        organism_ipr2ec = self.map_ipr_to_specific_ecs(organismName, organism_ipr2go)
                        organism_seqid2ec = self.combine_iprscan_raw_result_with_ipr2ec( organismName, organism_ipr2ec, ipr_raw_file)
                        
                       
                        if os.path.exists(ipr_raw_file) and os.path.exists(organism_seqid2ec):
    
    
                            organism_raw_final = NGS_Util.createFilePath(self.fungi_InterProScan_result, organismName + ".faa.raw.txt")                        
                            organism_IPR_final = NGS_Util.createFilePath(self.fungi_InterProScan_result, organismName + ".faa.IPR.final.txt")
                            
                            NGS_Util.copyFile(ipr_raw_file, organism_raw_final)
                            NGS_Util.copyFile(organism_seqid2ec, organism_IPR_final)


            orgListFile_fh.close() 
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""

