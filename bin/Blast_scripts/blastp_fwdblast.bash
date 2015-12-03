#!/bin/bash
IFS=$'\n'

#deine the input, output and software path
path_sequences=/users/jianhou/work/ncbi-blast-2.2.27+/sequences
path_software=/users/jianhou/work/ncbi-blast-2.2.27+/bin
path_makeblastdb=/users/jianhou/work/ncbi-blast-2.2.27+/results_makeblastdb
path_blast=/users/jianhou/work/ncbi-blast-2.2.27+/results_fwdblast

for name in $(ls $path_sequences -1)
do
  echo $name
  $path_software/blastp -db $path_software/uiprotdb -query $path_sequences/$name -outfmt 6 -out $path_blast/$name.2uniprot.txt -evalue 10
done




