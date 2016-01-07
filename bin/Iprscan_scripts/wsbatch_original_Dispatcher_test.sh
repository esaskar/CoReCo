#!/bin/bash
# ======================================================================
# Example script for submitting a set of single input per file jobs to 
# a Dispatcher based EMBL-EBI web service.
#
# See:
# http://www.ebi.ac.uk/Tools/webservices/
# http://www.ebi.ac.uk/Tools/webservices/help/faq
# http://www.ebi.ac.uk/Tools/webservices/tutorials/01_intro
# http://www.ebi.ac.uk/Tools/webservices/tutorials/07_workflows
# ----------------------------------------------------------------------
# NB: this script illustrates a possible submission workflow for 
# handling the submission of multiple concurrent jobs. In order to use 
# this script it will have to be modified to match your environment and 
# updated to run the required web service client with the required 
# options.
# ======================================================================
### Defaults ###
# User e-mail address. This *must* be a vaild e-mail address.
userEmail="jian102300@sina.com"
# Max number of simultaneous jobs.
declare -i MAX=1     # -i mean integer--------->declare OPTION(s) VARIABLE=value
# Executable for the web service client.
exec="python iprscan_suds.py"
#exec="java -Djava.ext.dirs=lib -jar bin/WSMaxsprout.jar"
# Command line options for the client.
# NB: the 'email' and 'async' options are required.
command_line="--email $userEmail";       # here $CHARACTER means get the valuse of this char

### Command-line ###
# You need to specify a directory to process
if [ $# -gt 0 ]; then          # $# means count the number after bashfile (.sh)
    todo="$1"                  # get the first argument behind .sh;if $0, get the .sh file name
else
    echo "Error: a directory to process must be specified" 1>&2
    exit 1
fi


### Run jobs ###
# Initialise tracking varibles.
declare -i pend=0      # Number of pending/running jobs.
# Remember where we are...
home=`pwd`
# Get list of files to process.
cd $todo
files=`ls`
# Move back to working directory.
cd $home

# For each file to process...
for seq in $files; do
    # Submit the job for the file.
    jobid=`${exec} ${command_line} ${todo}/${seq}`
   # echo "jobid------>$jobid"
    echo "[INFO] $jobid for $seq"
    # Record the job identifier.
   # PENDING="${PENDING} ${jobid}"
    # Update the pending/running job counter.
    pend=pend+1
    echo "pend is---->$pend"
  
done # for all files


