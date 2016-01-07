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


sys.path.append(ScriptsDir.GTGScripts)

class MetabolicReconstructionPipeline_GTG:

    nrdb40_blast_db  = ""
    
    ec_files          = ""     #need to be initialized     

    ngsBlast          = NGS_Blast.NGS_Blast()
    
    orgListFile       = ""  #Name of file containing Organisms List                  #need to be initialized     
    orgFastaDir       = ""  #Directory path containing organisms fasta sequences     #need to be initialized     

    orgBlastDBDir     = ""  #Directory path for organisms BLASTable databases        #need to be initialized     
    orgBlastDustDir   = ""  #Directory path for organisms DUST  files                #need to be initialized     

    
    orgGTGDatabaseDir   = ""  #Directory containing GTG input                        #need to be initialized     
    
    orgGTGBlastResDir   = ""  #Directory to store blast outputs                        #need to be initialized     
    GTGFungiBestHitsDir = ""  #Directory to store joint blast rectified outputs        #need to be initialized
    GTGKNNDir      = ""  #Directory to store joint blast rectified outputs        #need to be initialized
    
    CAA1Dir = ""     #need to be initialized     
    
    nids_up = ""     #need to be initialized
    
    seq_org_list = ""  #need to be absolute
    
    numberNearestHits = 50

    blastEValue       = 1e-10

    def initialize(self, nrdb40_blast_db, ec_files, orgListFile, orgFastaDir,  orgBlastDBDir, orgBlastDustDir, orgGTGBlastResDir, GTGBestHitsDir, GTGKNNDir, CAA1Dir, nids_up, seq_org_list, numberNearestHits, orgGTGDatabaseDir, blastEValue):
    
	    self.nrdb40_blast_db     = nrdb40_blast_db
	    
	    self.ec_files            = ec_files
    
	    self.orgListFile         = orgListFile
	    self.orgFastaDir         = orgFastaDir
	
	    self.orgBlastDBDir       = orgBlastDBDir
	    self.orgBlastDustDir     = orgBlastDustDir  
	
	    self.orgGTGDatabaseDir   = orgGTGDatabaseDir
	    self.orgGTGBlastResDir   = orgGTGBlastResDir
	    self.GTGBestHitsDir = GTGBestHitsDir
	    self.GTGKNNDir      = GTGKNNDir
	    
	    self.CAA1Dir             = CAA1Dir
	    self.nids_up             = nids_up
	    self.seq_org_list        = seq_org_list
	    self.numberNearestHits   = numberNearestHits
	    
	    self.blastEValue         = blastEValue
        
    def blast_org_vs_nr40_blast_formatted_11(self, organismName):
        print "blast_org_vs_nr40_blast_formatted_11: " + organismName
        org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName + ".faa")
	    org_vs_nr40BlastDB_f11 = NGS_Util.createFilePath(self.orgGTGBlastResDir, organismName + ".nrdb40_v2.txt")
	    self.ngsBlast.blastP(self.nrdb40_blast_db, org_fasta,  11, org_vs_nr40BlastDB_f11, self.blastEValue)
	    return org_vs_nr40BlastDB_f11

    def blast_org_vs_nr40_blast_formatted_6(self, organismName, org_vs_nr40BlastDB_f11):
        print "blast_org_vs_nr40_blast_formatted_6: " + organismName
        org_vs_nr40BlastDB_f6 = NGS_Util.createFilePath(self.orgGTGBlastResDir, organismName + ".nrdb40_v2_6.txt")
        self.ngsBlast.blastFormatter(org_vs_nr40BlastDB_f11, 6, org_vs_nr40BlastDB_f6)
        return org_vs_nr40BlastDB_f6

    def extract_seq_fmt11(self, organismName, org_vs_nr40BlastDB_f11): 	    #(2) extract query information from blast fmt11. : .part1
        print "extract_seq_fmt11: " + organismName

        org_vs_nr40BlastDB_f11_part1 = NGS_Util.createFilePath(self.orgGTGBlastResDir, organismName + ".nrdb40_v2.part1")
        call = "python " + ScriptsDir.GTGScripts_extract_seq_fmt11 + " " + org_vs_nr40BlastDB_f11 + " " + org_vs_nr40BlastDB_f11_part1
        NGS_Util.executeCall(call)
        return org_vs_nr40BlastDB_f11_part1

    def extract_start_len_fmt11(self, organismName, org_vs_nr40BlastDB_f11): #(3) extract start, len and subject name from fmt11 : .part2
        print "extract_start_len_fmt11: " + organismName
        org_vs_nr40BlastDB_f11_part2 = NGS_Util.createFilePath(self.orgGTGBlastResDir, organismName + ".nrdb40_v2.part2")
        call = "python " + ScriptsDir.GTGScripts_extract_start_len_fmt11 + " " + org_vs_nr40BlastDB_f11 + " " + org_vs_nr40BlastDB_f11_part2
        NGS_Util.executeCall(call)
        return org_vs_nr40BlastDB_f11_part2

    def extract_combine_seq_start_len_fmt11(self, organismName, org_vs_nr40BlastDB_f11_part1, org_vs_nr40BlastDB_f11_part2): #(4) combine the result from previous two steps
    
        print "extract_combine_seq_start_len_fmt11: " + organismName
        org_vs_nr40BlastDB_f11_part1_part2_result = NGS_Util.createFilePath(self.orgGTGBlastResDir, organismName + ".nrdb40_v2.part1.part2.result")
        call = "python " + ScriptsDir.GTGScripts_extract_combine_seq_start_len_fmt11 + " " + org_vs_nr40BlastDB_f11_part1 + " " + org_vs_nr40BlastDB_f11_part2 + " " + org_vs_nr40BlastDB_f11_part1_part2_result
        NGS_Util.executeCall(call)
        return org_vs_nr40BlastDB_f11_part1_part2_result

    def reform(self, organismName, org_vs_nr40BlastDB_f6, org_vs_nr40BlastDB_f11_part1_part2_result): #(5) extract and reform based on $name.nrdb40_v2_6.txt and $name.nrdb40_v2.txt: should give you the sample output below
    
        print "reform: " + organismName
        org_vs_nr40BlastDB_result_final = NGS_Util.createFilePath(self.orgGTGBlastResDir, organismName + ".nrdb40.result.final")
        call = "python " + ScriptsDir.GTGScripts_reform + " " + org_vs_nr40BlastDB_f6  + " " + org_vs_nr40BlastDB_f11_part1_part2_result + " " + org_vs_nr40BlastDB_result_final
        NGS_Util.executeCall(call)
        return org_vs_nr40BlastDB_result_final
	    
    def extract_best_hit(self, organismName, org_vs_nr40BlastDB_result_final): #(6) Extract the best hit for each query seq.
        print "extract_best_hit: " + organismName
        org_vs_nr40BlastDB_best_hit = NGS_Util.createFilePath(self.GTGBestHitsDir, organismName + ".nrdb40.best_hit")
        call = "python " + ScriptsDir.GTGScripts_extract_best_hit + " " + org_vs_nr40BlastDB_result_final  + " " + org_vs_nr40BlastDB_best_hit

        NGS_Util.executeCall(call)
        return org_vs_nr40BlastDB_best_hit

    def extract_gtg(self, organismName, org_vs_nr40BlastDB_best_hit):  #(7) Extract gtg feature for each org
        print "extract_gtg: " + organismName
        org_gtg = NGS_Util.createFilePath(self.GTGBestHitsDir, organismName + ".gtg")
        call = "cut -f 3,6-7,13-16 " + org_vs_nr40BlastDB_best_hit + " |time python " + ScriptsDir.GTGScripts_gtg_attributes_mod +  " 1 " + self.orgGTGDatabaseDir + " > " + org_gtg
        NGS_Util.executeCall(call)
        return org_gtg
 
    def gtgknn(self, organismName, org_gtg, numberNearestHits): #Find closest sequences in GTG space
        print "gtgknn: " + organismName
        org_gtg_knn = NGS_Util.createFilePath(self.GTGBestHitsDir, organismName + ".gtg.knn")
        call = "python " + ScriptsDir.GTGScripts_gtgknn + " " + org_gtg + " " + self.CAA1Dir + " " + self.nids_up + " " + str(numberNearestHits) + " " + org_gtg_knn
        NGS_Util.executeCall(call)
        return org_gtg_knn

    def reform_knn(self, organismName, org_gtg_knn): # (9) Add org and ecs
        print "reform_knn: " + organismName
        org_gtg_knn_final = NGS_Util.createFilePath(self.GTGKNNDir, organismName + ".gtg.knn")
        call = "python " + ScriptsDir.GTGScripts_reform_knn + " " +self.seq_org_list + " " + self.ec_files + " " + org_gtg_knn + " " + org_gtg_knn_final
        NGS_Util.executeCall(call)
        return org_gtg_knn_final
    
    def getGTGScore(self):
        orgListFile_fh = open(self.orgListFile)
        for line in orgListFile_fh:
            organismNameID, organismName = line.strip().split()
            org_gtg_knn_final = NGS_Util.createFilePath(self.GTGKNNDir, organismNameID + ".gtg.knn")
            if not os.path.exists(org_gtg_knn_final):

        print "getGTGScore : " +organismName

        org_vs_nr40BlastDB_f11 = self.blast_org_vs_nr40_blast_formatted_11( organismName)
        org_vs_nr40BlastDB_f6  = self.blast_org_vs_nr40_blast_formatted_6( organismName, org_vs_nr40BlastDB_f11)

        org_vs_nr40BlastDB_f11_part1 = self.extract_seq_fmt11( organismName, org_vs_nr40BlastDB_f11)
        org_vs_nr40BlastDB_f11_part2 = self.extract_start_len_fmt11( organismName, org_vs_nr40BlastDB_f11)
        org_vs_nr40BlastDB_f11_part1_part2_result = self.extract_combine_seq_start_len_fmt11( organismName, org_vs_nr40BlastDB_f11_part1, org_vs_nr40BlastDB_f11_part2)
        org_vs_nr40BlastDB_result_final = self.reform( organismName, org_vs_nr40BlastDB_f6, org_vs_nr40BlastDB_f11_part1_part2_result)
        org_vs_nr40BlastDB_best_hit     = self.extract_best_hit( organismNameID, org_vs_nr40BlastDB_result_final)
        org_gtg                         = self.extract_gtg( organismNameID, org_vs_nr40BlastDB_best_hit)
        org_gtg_knn                     = self.gtgknn( organismNameID, org_gtg, self.numberNearestHits)
        org_gtg_knn_final               = self.reform_knn( organismNameID, org_gtg_knn)
        orgListFile_fh.close() 

            
