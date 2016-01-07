#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""


import os
import sys
import traceback
sys.path.append("model_reconstruction_pipeline")
import NGS_Util

projectDir = "/project_path/coreco"

##########################################################################################     Blast Toolkit      ##########################################################################################

BlastDir     = "Tools/ncbi-blast-2.2.30+/bin/"
BlastDBDir   = "Tools/ncbi-blast-2.2.30+/db/"   #Needs to be set the user#
BlastDustDir = "Tools/ncbi-blast-2.2.30+/db/"   #Needs to be set the user#

##########################################################################################     Blast Toolkit      ##########################################################################################

IprscanDir   = "Tools/iprscan/bin/"

############################################################################################################################################################################################################

Java_PATH    = "/usr/java/latest/bin/"
LibSBML_PATH = "/share/apps/local/lib64/python2.6/site-packages/"

############################################################################################################################################################################################################

fastaSplitNDir  = "Path to the program to cut the fasta file in pieces called fastasplitN (Author: David Mathog, Biology Division, Caltech")
fastaSplitN     = NGS_Util.createFilePath(fastaSplitNDir,"fastasplitn")

###############################################################################################################################################################################################



projectBinDir  = NGS_Util.createDirectoryPath(projectDir,"bin")

BLASTScripts = NGS_Util.createDirectoryPath(projectDir,"Blast_scripts")
BlastScripts_blast=BLASTScripts + "blast.bash"
BlastScripts_blastp_bwdblast=BLASTScripts + "blastp_bwdblast.bash"
BlastScripts_blastp_fwdblast=BLASTScripts + "blastp_fwdblast.bash"
BlastScripts_makeblastdb_AUTO=BLASTScripts + "makeblastdb_AUTO.bash"
BlastScripts_makedustfile_AUTO=BLASTScripts + "makedustfile_AUTO.bash"
BlastScripts_buildBlastResult=BLASTScripts + "buildBlastResult.py"
BlastScripts_combineAllBlasts=BLASTScripts + "combineAllBlasts.py"
BlastScripts_getEcs=BLASTScripts + "getEcs.py"
BlastScripts_rectify_blastresult=BLASTScripts + "rectify_blastresult.py"


GTGScripts = NGS_Util.createDirectoryPath(projectDir,"GTG_scripts")
GTGScripts_reform_knn= GTGScripts + "reform_knn.py"
GTGScripts_reform= GTGScripts + "reform.py"
GTGScripts_linebuffer= GTGScripts + "linebuffer.py"
GTGScripts_gtgknn= GTGScripts + "gtgknn.py"
GTGScripts_gtg_attributes_mod= GTGScripts + "gtg_attributes_mod.py"
GTGScripts_extract_start_len_fmt11= GTGScripts + "extract_start_len_fmt11.py"
GTGScripts_extract_seq_fmt11= GTGScripts + "extract_seq_fmt11.py"
GTGScripts_extract_nids_from_uniprot= GTGScripts + "extract_nids_from_uniprot.py"
GTGScripts_extract_combine_seq_start_len_fmt11= GTGScripts + "extract_combine_seq_start_len_fmt11.py"
GTGScripts_extract_best_hit= GTGScripts + "extract_best_hit.py"
GTGScripts_buildGTGindex= GTGScripts + "buildGTGindex.py"

IPRScanScripts = NGS_Util.createDirectoryPath(projectDir,"Iprscan_scripts")
IPRScanScripts_wsbatch_original_Dispatcher_test = IPRScanScripts + "wsbatch_original_Dispatcher_test.sh"
IPRScanScripts_iprscan_suds = IPRScanScripts + "iprscan_suds.py"
IPRScanScripts_ipr2go = IPRScanScripts + "ipr2go.py"
IPRScanScripts_ipr_reform_ecs = IPRScanScripts + "ipr_reform_ecs.bash"
IPRScanScripts_ipr_bash = IPRScanScripts + "ipr_bash.bash"
IPRScanScripts_get_seq_org_list = IPRScanScripts + "get_seq_org_list.py"
IPRScanScripts_get_interpro_ecs = IPRScanScripts + "get_interpro.ecs.py"
IPRScanScripts_fsplit = IPRScanScripts + "fsplit.sh"
IPRScanScripts_combineIPRwithECs = IPRScanScripts + "combineIPRwithECs.py"
IPRScanScripts_combine_xml_out = IPRScanScripts + "combine_xml_out.bash"

ModelTrainingScripts = NGS_Util.createDirectoryPath(projectDir,"model_training_scripts")
ModelTrainingScripts_visualize_ec_scores = ModelTrainingScripts  + "visualize_ec_scores.py"
ModelTrainingScripts_tree = ModelTrainingScripts  + "tree.py"
ModelTrainingScripts_run_job = ModelTrainingScripts  + "run_job.sh"
ModelTrainingScripts_roc = ModelTrainingScripts  + "roc.py"
ModelTrainingScripts_plot_reaction_scores = ModelTrainingScripts  + "plot_reaction_scores.py"
ModelTrainingScripts_plot = ModelTrainingScripts  + "plot.py"
ModelTrainingScripts_merge_scores = ModelTrainingScripts  + "merge_scores.py"
ModelTrainingScripts_import_data = ModelTrainingScripts  + "import_data.py"
ModelTrainingScripts_fitch = ModelTrainingScripts  + "fitch.py"
ModelTrainingScripts_extract_ecs_from_iprscan = ModelTrainingScripts  + "extract_ecs_from_iprscan.py"
ModelTrainingScripts_estimate_mutation_probability = ModelTrainingScripts  + "estimate_mutation_probability.py"
ModelTrainingScripts_estimate_cpds_wrapper = ModelTrainingScripts  + "estimate_cpds_wrapper.py"
ModelTrainingScripts_estimate_cpds = ModelTrainingScripts  + "estimate_cpds.py"
ModelTrainingScripts_digraph = ModelTrainingScripts  + "digraph.py"
ModelTrainingScripts_computeMergedScores = ModelTrainingScripts  + "computeMergedScores.py"
ModelTrainingScripts_computeECscores = ModelTrainingScripts  + "computeECscores.py"
ModelTrainingScripts_computeBlastPvalues = ModelTrainingScripts  + "computeBlastPvalues.py"
ModelTrainingScripts_compute_reaction_scores = ModelTrainingScripts  + "compute_reaction_scores.py"
ModelTrainingScripts_combined_density = ModelTrainingScripts  + "combined_density.py"
ModelTrainingScripts_bayesnet = ModelTrainingScripts  + "bayesnet.py"


ReconstructionScripts = NGS_Util.createDirectoryPath(projectDir,"reconstruction_scripts")
ReconstructionScripts_reco_dir =  ReconstructionScripts + "reco-dir"
ReconstructionScripts_reco_dir_postprocess =  ReconstructionScripts + "reco-dir-postprocess.sh"
ReconstructionScripts_reco_dir_cluster =  ReconstructionScripts + "reco-dir_cluster.sh"
ReconstructionScripts_reco_dir_postprocess_cluster =  ReconstructionScripts + "reco-dir-postprocess_cluster.sh"

keggParsingScripts =  NGS_Util.createDirectoryPath(projectDir,"kegg-parsing")
keggParsingScripts_build_kegg_no_general = keggParsingScripts + "build_kegg_no_general.sh"

modelReconstructionPipelineScripts = NGS_Util.createDirectoryPath(projectBinDir,"model_reconstruction_pipeline")
ClusterBlast    = modelReconstructionPipelineScripts  + "ClusterBlast.sh"
ClusterIprscan  = modelReconstructionPipelineScripts  + "ClusterIprscan.sh"
ClusterGTGBlast = modelReconstructionPipelineScripts  + "ClusterGTGBlast.sh"
ClusterPipeline = modelReconstructionPipelineScripts  + "ClusterPipeline.sh"
