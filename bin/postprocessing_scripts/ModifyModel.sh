#$ -S /bin/bash
#$ -N ModifyModel
#$ -j y
#$ -cwd
#$ -pe smp 4
#$ -R y

# Post processing step independent on the reconstruction pipeline. 
# Once the models are created with CoReCo pipeline, this scripts allows you to add reactions,
# the biomass equation and new bounds depending on the choosen parameters. Example imput files can be found in the folder "example_input_files". 

path_models=example_input_files/models
path_biomass=example_input_files/biomass.csv
path_bounds=example_input_files/harvestedbounds.csv
path_reactions=example_input_files/reactions2add.csv
output_path=example_input_files/output  

for name in $(ls $path_models -1)
do
  echo $name
  python ModifyModel.py -m  $path_models/$name -r $path_reactions  -b $path_biomass -k $path_bounds -o $output_path/$name
done


