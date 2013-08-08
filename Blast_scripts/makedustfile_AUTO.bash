#! /bin/bash
IFS=$'\n'

#deine the input, output and software path
path_sequences=/users/jianhou/work/ncbi-blast-2.2.27+/sequences
path_software=/users/jianhou/work/ncbi-blast-2.2.27+/bin
path_makeblastdb=/users/jianhou/work/ncbi-blast-2.2.27+/results_makeblastdb


for name in $(ls $path_sequences -1)
do
  echo $name
  $path_software/segmasker -in $path_sequences/$name -infmt fasta -outfmt maskinfo_asn1_bin -parse_seqids -out $path_makeblastdb/$name.asnb
done


