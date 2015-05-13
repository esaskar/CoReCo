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

sys.path.append(ScriptsDir.BLASTScripts)

from buildBlastResult  import combineBlasts

class MetabolicReconstructionPipeline_Blast:

    uniprot_fasta     = ""     #need to be initialized     
    uniprot_dust      = ""
    uniprot_blast_db  = ""

    ec_files          = ""     #need to be initialized     

    uniprot_sprot_dat = ""     #need to be initialized      
    
    ngsBlast          = NGS_Blast.NGS_Blast()
    
    orgListFile       = ""  #Name of file containing Organisms List                  #need to be initialized     
    orgFastaDir       = ""  #Directory path containing organisms fasta sequences     #need to be initialized     

    orgBlastDBDir     = ""  #Directory path for organisms BLASTable databases        #need to be initialized     
    orgBlastDustDir   = ""  #Directory path for organisms DUST  files                #need to be initialized     

    orgBlastResDir    = ""  #Directory to store blast outputs                        #need to be initialized     
    jointBlastDir     = ""  #Directory to store joint blast rectified outputs        #need to be initialized     

    def initialize(self, uniprot_fasta, ec_files, uniprot_sprot_dat, uniprot_dust, uniprot_blast_db, orgListFile, orgFastaDir , orgBlastDBDir, orgBlastDustDir, orgBlastResDir, jointBlastDir):
    
        self.uniprot_fasta     = uniprot_fasta
        self.uniprot_dust      = uniprot_dust
        self.uniprot_blast_db  = uniprot_blast_db

        self.ec_files          = ec_files

        self.uniprot_sprot_dat = uniprot_sprot_dat


        self.orgListFile       = orgListFile
        self.orgFastaDir       = orgFastaDir

        self.orgBlastDBDir     = orgBlastDBDir
        self.orgBlastDustDir   = orgBlastDustDir

        self.orgBlastResDir    = orgBlastResDir
        self.jointBlastDir     = jointBlastDir
        
    def makeBlastDB(self, organismName):
        print "Make Blast Database: " + organismName
        org_fasta    = NGS_Util.createFilePath(self.orgFastaDir, organismName+".faa")
        org_dust     = NGS_Util.createFilePath(self.orgBlastDustDir, organismName+"_dust.asnb")
        org_blast_db = NGS_Util.createFilePath(self.orgBlastDBDir, organismName)
        if os.path.exists(org_fasta):
            if not os.path.exists(org_blast_db + ".phd") and not os.path.exists(org_blast_db + ".psq"):
                self.ngsBlast.makeProteinBlastDBFromDustFile(org_fasta,org_dust,org_blast_db)
            return org_blast_db

    def blast_org_vs_uniprot(self, organismName):
        print "org_vs_uniprot_blast: " + organismName
        org_fasta = NGS_Util.createFilePath(self.orgFastaDir, organismName+".faa")
        org_vs_UniprotBlastDB =  NGS_Util.createFilePath(self.orgBlastResDir, organismName+"-vs-up.blast")
        self.ngsBlast.blastP(self.uniprot_blast_db,org_fasta, 6 , org_vs_UniprotBlastDB, 10)
        return org_vs_UniprotBlastDB


    def blast_uniprot_vs_org(self, organismName):
        print "uniprot_vs_org_blast Blast: " + organismName
        org_blast_db = NGS_Util.createFilePath(self.orgBlastDBDir, organismName)
        Uniprot_vs_orgBlastDB = NGS_Util.createFilePath(self.orgBlastResDir, "up-vs-" + organismName+ ".blast")
        self.ngsBlast.blastP(org_blast_db, self.uniprot_fasta,6, Uniprot_vs_orgBlastDB, 10)
        return Uniprot_vs_orgBlastDB

    def combineBlast(self, organismName, org_vs_UniprotBlastDB, Uniprot_vs_orgBlastDB):
        print "Combine Blast: " + organismName

        orgJointBlast            = NGS_Util.createFilePath(self.orgBlastResDir, organismName + ".joint.blast")
        org_vs_UniprotBlastDB_fh = open(org_vs_UniprotBlastDB)
        Uniprot_vs_orgBlastDB_fh = open(Uniprot_vs_orgBlastDB)
        ec_files_fh              = open(self.ec_files)
        orgJointBlast_fh         = open(orgJointBlast, "w")

        combineBlasts(ec_files_fh, org_vs_UniprotBlastDB_fh, Uniprot_vs_orgBlastDB_fh, orgJointBlast_fh)
 
        org_vs_UniprotBlastDB_fh.close()
        Uniprot_vs_orgBlastDB_fh.close()
        ec_files_fh.close()
        orgJointBlast_fh.close()

        return orgJointBlast       

    def rectifyBlast(self, organismName, orgJointBlast):
        print "Rectify Blast: " + organismName

        orgRectifyBlast = NGS_Util.createFilePath(self.jointBlastDir, organismName + ".joint.blast")
        call = "python " + ScriptsDir.BlastScripts_rectify_blastresult +  " " + orgJointBlast + " " + orgRectifyBlast
        NGS_Util.executeCall(call)
        return orgRectifyBlast

    def getBlastScore(self):
        orgListFile_fh = open(self.orgListFile)

        for line in orgListFile_fh:
            if line.startswith("#"):
                continue

            organismNameID, organismName = line.strip().split()
            orgRectifyBlast = NGS_Util.createFilePath(self.jointBlastDir, organismName + ".joint.blast")

            if not os.path.exists(orgRectifyBlast):
                print "getBlastScore:" + organismName
                org_blast_db = self.makeBlastDB(organismName)
                org_vs_UniprotBlastDB = self.blast_org_vs_uniprot(organismName)
                Uniprot_vs_orgBlastDB = self.blast_uniprot_vs_org(organismName)
                if (org_vs_UniprotBlastDB != "" and Uniprot_vs_orgBlastDB != ""):
                    orgJointBlast = self.combineBlast(organismName, org_vs_UniprotBlastDB, Uniprot_vs_orgBlastDB)
                    if (orgJointBlast != ""):
                        orgRectifyBlast = self.rectifyBlast(organismName, orgJointBlast)

        orgListFile_fh.close() 
     


