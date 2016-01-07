#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""

import os
import sys
import subprocess
import operator
import traceback
import time
import threading
import glob

import NGS_Util
import NGS_IPRScan

import MetabolicReconstructionPipeline_SplitFasta

sys.path.append("..")
import ScriptsDir

sys.path.append(ScriptsDir.IprscanDir)

class MetabolicReconstructionPipeline_IPRScan:

    ngsIPRScan       = NGS_IPRScan.NGS_IPRScan()

    orgIPRScanDir             = ""
    fungi_InterProScan_result = ""

    ec2go            = ""  #need to be initialized     
    
    orgListFile      = ""  #Name of file containing Organisms List                  #need to be initialized     
    orgFastaDir      = ""  #Directory path containing organisms fasta sequences     #need to be initialized     

    seq_org_list     = ""  #need to be absolute
    
    numberOfFragments = 10

    splitFasta = MetabolicReconstructionPipeline_SplitFasta.MetabolicReconstructionPipeline_SplitFasta()

    def initialize(self, orgListFile, orgFastaDir,  orgIPRScanDir, fungi_InterProScan_result, ec2go, seq_org_list, numberOfFragments): 
        self.orgListFile         = orgListFile
        self.orgFastaDir         = orgFastaDir
        self.orgIPRScanDir             = orgIPRScanDir
        self.fungi_InterProScan_result = fungi_InterProScan_result
        self.ec2go            = ec2go
        self.seq_org_list     = seq_org_list
	    self.numberOfFragments = numberOfFragments
   
    def splitFiles(self, organismName):
        org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName+".faa")
        self.splitFasta.splitOrganismDataFile(organismName, org_fasta, self.numberOfFragments)
        return self.splitFasta.organismSplitDataDir

    def raw_split_IPRScan(self, organismName, organismSplitFile, splitNameIndex):
        print "raw_split_IPRScan: " + organismName + " " + organismSplitFile
	    ipr_raw_file_split = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + "_split_" +  str(splitNameIndex) + ".ipr.raw")
        if not os.path.exists(ipr_raw_file_split):
            self.ngsIPRScan.protein_iprscan_to_raw_output(organismSplitFile, ipr_raw_file_split)
        return ipr_raw_file

    def rawIPRScan(self, organismName, org_ipr_split_dir):
        print "rawIPRScan: " + organismName
        for fragment in range(self.numberOfFragments):
            org_ipr_split_file = NGS_Util.createFilePath(self.splitFasta.organismSplitDataDir,organismName + "_" +  str(fragment+1) )
            self.raw_split_IPRScan(organismName, org_ipr_split_file, str(fragment+1))
        ipr_raw_file_split = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + "_split_*")            
        ipr_raw_file = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + ".ipr.raw")
        call = "cat " + ipr_raw_file_split + " > " + ipr_raw_file
	    NGS_Util.executeCall(call)
        return ipr_raw_file
    
    def rawIPRScanToXMlOutput(self, organismName, ipr_raw_file):
        print "rawIPRScanToXMlOutput: " + organismName
	    ipr_xml_file = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + ".xml")
	    self.ngsIPRScan.convert_raw_xml(ipr_raw_file, ipr_xml_file)
	    return ipr_xml_file
        
    def getIPRScanScore(self):
        orgListFile_fh = open(self.orgListFile)

        for line in orgListFile_fh:
            organismNameID, organismName = line.strip().split()
            organism_IPR_final = NGS_Util.createFilePath(self.fungi_InterProScan_result, organismName + ".faa.IPR.final.txt")

            if not os.path.exists(organism_IPR_final):
                print "getIPRScanScore : " + organismName
                org_ipr_split_dir = self.splitFiles(organismName)
                ipr_raw_file = self.rawIPRScan(organismName,org_ipr_split_dir)
                ipr_xml_file = self.rawIPRScanToXMlOutput( organismName, ipr_raw_file)

                organism_ipr2go = self.extract_ipr2go_based_on_xml(organismName, ipr_xml_file)
                organism_ipr2ec = self.map_ipr_to_specific_ecs(organismName, organism_ipr2go)
                organism_seqid2ec = self.combine_iprscan_raw_result_with_ipr2ec( organismName, organism_ipr2ec, ipr_raw_file)

                if os.path.exists(ipr_raw_file) and os.path.exists(organism_seqid2ec):

                    organism_raw_final = NGS_Util.createFilePath(self.fungi_InterProScan_result, organismName + ".faa.raw.txt")                        
                    organism_IPR_final = NGS_Util.createFilePath(self.fungi_InterProScan_result, organismName + ".faa.IPR.final.txt")

                    NGS_Util.copyFile(ipr_raw_file, organism_raw_final)
                    NGS_Util.copyFile(organism_seqid2ec, organism_IPR_final)

        orgListFile_fh.close() 

    def extract_ipr2go_based_on_xml(self,organismName, ipr_xml_file):
        print "extract_ipr2go_based_on_xml: " + organismName
	    organism_ipr2go = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + "_ipr2go.txt")
	    call = "python " + ScriptsDir.IPRScanScripts_ipr2go + " " + ipr_xml_file + " " + organism_ipr2go
	    NGS_Util.executeCall(call)
	    return organism_ipr2go

    def map_ipr_to_specific_ecs(self,organismName, organism_ipr2go):
        print "map_ipr_to_specific_ecs: " + organismName
	    organism_ipr2ec = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + "_ipr2ec.txt")
	    call = "python " + ScriptsDir.IPRScanScripts_get_interpro_ecs + " " + self.ec2go + " " + organism_ipr2go + " " + organism_ipr2ec
	    NGS_Util.executeCall(call)
	    return organism_ipr2ec

    def combine_iprscan_raw_result_with_ipr2ec(self, organismName, organism_ipr2ec, ipr_raw_file):  ### to be changes new_seq_org_list -> seq_org_list
        print "combine_iprscan_raw_result_with_ipr2ec: " + organismName
	    organism_seqid2ec = NGS_Util.createFilePath(self.orgIPRScanDir, organismName + "_seqid2ec.txt")
	    call = "python " + ScriptsDir.IPRScanScripts_combineIPRwithECs + " " + organism_ipr2ec + " " + ipr_raw_file + " " + self.seq_org_list + " " + organism_seqid2ec
	    NGS_Util.executeCall(call)
	    return organism_seqid2ec
        