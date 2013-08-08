#!/bin/bash
IFS=$'\n'

#deine the input seq (the same as output dir) and software path
path_GTG=/users/jianhou/work/iprscan/50fungi         #Ashbya_gossypii.faa-python-suds-0.4	
path_software=/triton/ics/work/houj1/scripts
path_sequences=/users/jianhou/work/ncbi-blast-2.2.27+/sequences    #Ashbya_gossypii.faa
path_result=/users/jianhou/work/iprscan/50fungi_final	

for name in $(ls $path_sequences -1)
do
  echo $name
  cat $path_GTG/$name-python-suds-0.4/*.out.txt > $path_result/$name.raw.txt
  cat $path_GTG/$name-python-suds-0.4/*.xml > $path_result/$name.raw.xml 
done

wait

