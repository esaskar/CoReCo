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



class MetabolicReconstructionPipeline_PreProcess:

    uniprot_fasta     = "" #need to be initialized     
    uniprot_dust      = ""
    uniprot_blast_db  = ""

    nrdb40_fasta     = ""  #need to be initialized     
    nrdb40_dust      = ""
    nrdb40_blast_db  = ""


    ec_files          = ""  #need to be initialized     
    uniprot_sprot_dat = ""  #need to be initialized          
    
    
    orgListFile       = ""  #Name of file containing Organisms List                  #need to be initialized
    orgFastaDir       = ""  #Directory path containing organisms fasta sequences     #need to be initialized         
    seq_org_list      = ""  #Name of file containing Organisms List with GeneIdentifiers List #need to be initialized
    
    
    orgBlastDBDir     = ""  #Directory path for organisms BLASTable databases        #need to be initialized     
    orgBlastDustDir   = ""  #Directory path for organisms DUST  files                #need to be initialized     


    keggAtomMapsDataDir = ""

    ngsBlast          = NGS_Blast.NGS_Blast()
    

    def initialize(self, uniprot_fasta, uniprot_dust, uniprot_blast_db,  nrdb40_fasta, nrdb40_dust, nrdb40_blast_db, orgListFile, orgFastaDir, seq_org_list, orgBlastDBDir, orgBlastDustDir, ec_files, uniprot_sprot_dat, keggAtomMapsDataDir ):
    
        try:
        
            self.uniprot_fasta     = uniprot_fasta
	    self.uniprot_dust      = uniprot_dust
	    self.uniprot_blast_db  = uniprot_blast_db
	    
	    self.nrdb40_fasta      = nrdb40_fasta
	    self.nrdb40_dust       = nrdb40_dust
	    self.nrdb40_blast_db   = nrdb40_blast_db

            
            self.orgListFile       = orgListFile
	    self.orgFastaDir       = orgFastaDir   
            self.seq_org_list      = seq_org_list
        
            self.orgBlastDBDir     = orgBlastDBDir
            self.orgBlastDustDir   = orgBlastDustDir

            self.ec_files          = ec_files
            self.uniprot_sprot_dat = uniprot_sprot_dat

	    self.keggAtomMapsDataDir = keggAtomMapsDataDir

        except Exception:
            print traceback.print_exc()
    
    

    def makeUniprotBlastDB(self):
    
        try:
        
            print "makeUniprotBlastDB"
            

            if os.path.exists(self.uniprot_fasta):
            
                if not os.path.exists(self.uniprot_blast_db + ".phd") and not os.path.exists(self.uniprot_blast_db + ".psq"):
                
                    self.ngsBlast.makeProteinBlastDBFromDustFile(self.uniprot_fasta,self.uniprot_dust,self.uniprot_blast_db)

    
        except Exception:
            print traceback.print_exc()
    


    def preprocess_EC_UniprotSprot_data_files(self):
    
        try:
        
            print "preprocess_EC_UniprotSprot_data_files"

            if not os.path.exists( self.ec_files ) and os.path.exists(self.uniprot_sprot_dat):
    
                call = "python " + ScriptsDir.BlastScripts_getEcs + " " + self.uniprot_sprot_dat+ " " + self.ec_files +": AC, ID, EC"
            
                NGS_Util.executeCall(call)
        
    
        except Exception:
            print traceback.print_exc()
    
    
    def makeNrdb40BlastDB(self):
    
        try:
        
            print "makeNrdb40BlastDB"
            

            if os.path.exists(self.nrdb40_fasta):
            
                if not os.path.exists(self.nrdb40_blast_db + ".phd") and not os.path.exists(self.nrdb40_blast_db + ".psq"):

		    self.ngsBlast.makeProteinBlastDBFromDustFile(self.nrdb40_fasta, self.nrdb40_dust, self.nrdb40_blast_db)

        except Exception:

	    print traceback.print_exc()


       
    def create_new_seq_org_list(self):
    
        try:

            print "create_new_seq_org_list"
        
            orgListFile_fh = open(self.orgListFile)

            for orgLine in orgListFile_fh:
                
                organismID, organismName = orgLine.strip().split()
                
		seqOrgListFile_fh = open(self.seq_org_list)
		
		found = False
		
		for line in seqOrgListFile_fh:
			
			if "\t"+organismID in line:
				found =  True
				break
		
		seqOrgListFile_fh.close
		
		
		if not found:
		    
		    print "create_new_seq_org_list: " + organismName		    

		    org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName + ".faa")
		    
		    org_fasta_fh = open(org_fasta)

		    seqOrgListFile_fh = open(self.seq_org_list,"a")      #output file

		    for line in org_fasta_fh:

			if line.startswith(">"):

				seq_id = line.split(">")[1]

				if "|" in seq_id:
					id = seq_id.split("|")[0].strip()
				else:
					id = seq_id.strip()


#				if "|" in line:
#					id = line.split()[0].split("|")[1]
#				else:
#					id = line.split(" ")[0]
					
				seqOrgListFile_fh.write( id + "\t" +  organismID + "\n" )

                    seqOrgListFile_fh.close()
		    org_fasta_fh.close()


            orgListFile_fh.close() 
     
        except Exception:
            
            print traceback.print_exc()
            
        return ""

    
    
    def generateKeggData(self):
    
        try:
        
            print "generateKeggData"
            
            if not os.path.exists(self.keggAtomMapsDataDir):	    

		call = "bash " + ScriptsDir.keggParsingScripts_build_kegg_no_general
            
		NGS_Util.executeCall(call)
        
    

        except Exception:

	    print traceback.print_exc()




    def preProcess(self):
    
        try:
        
            self.makeUniprotBlastDB()
	    self.preprocess_EC_UniprotSprot_data_files()
            
            self.makeNrdb40BlastDB()
	    self.create_new_seq_org_list()

	    self.generateKeggData()

        except Exception:
            print traceback.print_exc()
        
