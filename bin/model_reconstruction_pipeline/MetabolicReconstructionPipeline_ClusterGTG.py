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
import MetabolicReconstructionPipeline_GTG

sys.path.append("..")
import ScriptsDir


class MetabolicReconstructionPipeline_ClusterGTG(MetabolicReconstructionPipeline_GTG.MetabolicReconstructionPipeline_GTG):
    
    def blast_org_vs_nr40_blast_formatted_11(self, organismName):
    
        try:
           
            org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName+".faa")

            clusterArrayCall =  "qsub -t 1-1 " + ScriptsDir.ClusterGTGBlast
            
            blastP = NGS_Util.createFilePath(ScriptsDir.BlastDir,"blastp")
            
            outfmt = str(11)
            
            org_vs_nr40BlastDB_f11 = NGS_Util.createFilePath(self.orgGTGBlastResDir, organismName + ".nrdb40_v2.txt")    

            call = clusterArrayCall + " "  + blastP  + " "  + self.nrdb40_blast_db + " "  + org_fasta + " " + outfmt + " " + org_vs_nr40BlastDB_f11  + " " + str(self.blastEValue)
                
            NGS_Util.executeCall(call)


            return org_vs_nr40BlastDB_f11
            	    
    
        except Exception:
            
            print traceback.print_exc()
            
        return ""


    def getGTGScore(self, mode):
    
        try:

	    counter = 0

            print "getGTGScore"
        
            orgListFile_fh = open(self.orgListFile)
	    
            for line in orgListFile_fh:
                
                organismNameID, organismName = line.strip().split()

                org_gtg_knn_final = NGS_Util.createFilePath(self.GTGKNNDir, organismNameID + ".gtg.knn")
		
                if not os.path.exists(org_gtg_knn_final):

		    print "getGTGScore : " + organismName
		    
		    if (mode == 1):
		    
			org_vs_nr40BlastDB_f11 = self.blast_org_vs_nr40_blast_formatted_11( organismName)
			
			counter += 1
			if (counter == 10):
			    time.sleep(5600) #wait for 30 minutes
			    counter = 0
			

		    elif(mode == 2):

	    		org_gtg_knn = NGS_Util.createFilePath(self.GTGBestHitsDir, organismNameID + ".gtg.knn")

			if os.path.exists(org_gtg_knn):
			    org_gtg_knn_final = self.reform_knn( organismNameID, org_gtg_knn)
			    
			else:

                            org_vs_nr40BlastDB_f11_part1_part2_result = NGS_Util.createFilePath(self.orgGTGBlastResDir, organismName + ".nrdb40_v2.part1.part2.result")
                            if not os.path.exists(org_vs_nr40BlastDB_f11_part1_part2_result):
                                org_vs_nr40BlastDB_f11 = NGS_Util.createFilePath(self.orgGTGBlastResDir, organismName + ".nrdb40_v2.txt")
                                org_vs_nr40BlastDB_f6  = self.blast_org_vs_nr40_blast_formatted_6( organismName, org_vs_nr40BlastDB_f11)
                                org_vs_nr40BlastDB_f11_part1 = self.extract_seq_fmt11( organismName, org_vs_nr40BlastDB_f11)
                                org_vs_nr40BlastDB_f11_part2 = self.extract_start_len_fmt11( organismName, org_vs_nr40BlastDB_f11)
                                org_vs_nr40BlastDB_f11_part1_part2_result = self.extract_combine_seq_start_len_fmt11( organismName, org_vs_nr40BlastDB_f11_part1, org_vs_nr40BlastDB_f11_part2)
                            else: 
                                print("\t %s already exists! Not recomputing" % (org_vs_nr40BlastDB_f11_part1_part2_result))
                            
                            org_vs_nr40BlastDB_result_final = NGS_Util.createFilePath(self.orgGTGBlastResDir, organismName + ".nrdb40.result.final")	    
                            if not os.path.exists(org_vs_nr40BlastDB_result_final):
                                org_vs_nr40BlastDB_result_final = self.reform( organismName, org_vs_nr40BlastDB_f6, org_vs_nr40BlastDB_f11_part1_part2_result)
                            else: 
                                print("\t %s already exists! Not recomputing" % (org_vs_nr40BlastDB_result_final))
			    
                            org_gtg = NGS_Util.createFilePath(self.GTGBestHitsDir, organismName + ".gtg")
			    if not os.path.exists(org_gtg):
                                org_vs_nr40BlastDB_best_hit     = self.extract_best_hit( organismNameID, org_vs_nr40BlastDB_result_final)
                                org_gtg                         = self.extract_gtg( organismNameID, org_vs_nr40BlastDB_best_hit)
                            else: 
                                print("\t %s already exists! Not recomputing" % (org_gtg))

                            org_gtg_knn = NGS_Util.createFilePath(self.GTGBestHitsDir, organismName + ".gtg.knn")
			    if not os.path.exists(org_gtg_knn):
                                org_gtg_knn                     = self.gtgknn( organismNameID, org_gtg, self.numberNearestHits)
                            else: 
                                print("\t %s already exists! Not recomputing" % (org_gtg_knn))
                                
                            org_gtg_knn_final = NGS_Util.createFilePath(self.GTGKNNDir, organismName + ".gtg.knn")
			    if not os.path.exists(org_gtg_knn_final):
                                org_gtg_knn_final               = self.reform_knn( organismNameID, org_gtg_knn)
                            else: 
                                print("\t %s already exists! Not recomputing" % (org_gtg_knn_final))
			  
            orgListFile_fh.close() 
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""


