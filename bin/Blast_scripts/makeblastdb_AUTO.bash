#! /bin/bash
IFS=$'\n'

#deine the input, output and software path
path_sequences=/users/jianhou/work/ncbi-blast-2.2.27+/sequences
path_software=/users/jianhou/work/ncbi-blast-2.2.27+/bin
path_makeblastdb=/users/jianhou/work/ncbi-blast-2.2.27+/results_makeblastdb


for name in $(ls $path_sequences -1)
do 

  $path_software/makeblastdb -in $path_sequences/$name -input_type fasta -dbtype prot -mask_data $path_makeblastdb/$name.asnb -parse_seqids -out $path_makeblastdb/$name.db

done






