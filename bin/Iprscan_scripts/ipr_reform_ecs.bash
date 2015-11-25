#!/bin/bash
IFS=$'\n'
#deine the input seq (the same as output dir) and software path
path_GTG=/users/jianhou/work/iprscan/50fungi_final        #Ashbya_gossypii.faa.raw.txt	Ashbya_gossypii.faa.raw.xml
path_software=/users/jianhou/work/iprscan/
path_sequences=/users/jianhou/work/ncbi-blast-2.2.27+/sequences    #Ashbya_gossypii.faa
path_result=/users/jianhou/work/iprscan/50fungi_final	

echo "extract ipr2go based on .xml ipr result"
for name in $(ls $path_sequences -1)
do
  echo $name
  python $path_software/ipr2go.py $path_GTG/$name.raw.xml $path_result/$name.2go.txt
done

wait

echo "map ipr to specific ecs based on ipr2go and ec2go"
for name in $(ls $path_sequences -1)
do
  echo $name
  python $path_software/get_interpro.ecs.py ec2go.txt $path_result/$name.2go.txt $path_result/$name.2ec.txt
done

wait

echo "combine iprscan_raw_result with ipr2ec.txt"
for name in $(ls $path_sequences -1)
do
  echo $name
  python $path_software/combineIPRwithECs.py $path_result/$name.2ec.txt $path_result/$name.raw.txt $path_software/seq_org_list.txt $path_result/$name.IPR.final.txt
done

wait

echo "Finish!"
