#!/bin/bash

#path explain
#sequences:50fungi .fasta sequence
#split: split .fasta into single .fasta
#scripts: save the running script .bash file
#result: save the 50fungi iprscan results
#suds: python-suds-0.4 path
path_scripts=/users/jianhou/work/iprscan/ipr_scripts
path_split=/users/jianhou/work/iprscan/split_sequences
path_sequences=/users/jianhou/work/ncbi-blast-2.2.27+/sequences
path_result=/users/jianhou/work/iprscan/50fungi
path_suds=/home/jianhou/Desktop/work/python-suds-0.4

# make dirctory of each species'name under path_split
# split the original .faa sequences into single and save it into $name-split
for name in $(ls $path_sequences -1)
do
  echo "----> $name-split"
  bash $path_scripts/fsplit.sh $path_sequences/$name $path_split/$name-split 
  echo "$name done"
done

# make 50 fungi python-suds-0.4
# copy all suds files in each fungi species 
for name in $(ls $path_sequences -1)
do 
  cp -r $path_suds $path_result/$name-python-suds-0.4
done

