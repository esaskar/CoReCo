#!/bin/bash
# 
# $1  data dir 
# $2  param-accept
# $3  param-reject
# $4  species
# $5  label-outputdir-with-params
# $8  bound file
# $9  ec2rxn file
# ${10}  exchange file

if [ -z "$1" ]; then
   echo "Specify data dir";
   exit;
fi

if [ -z "$2" ]; then
   echo "Specify accept threshold";
   exit;
fi

if [ -z "$3" ]; then
   echo "Specify reject threshold";
   exit;
fi

if [ -z "$4" ]; then
   echo "Specify species";
   exit;
fi


if [ -z "$5" ]; then
   echo "Specify Reconstruction Scripts Dir";
   exit;
fi


if [ -z "$6" ]; then
   echo "Specify Kegg Dir";
   exit;
fi

if [ -z "$7" ]; then
   echo "Specify taxonomy file";
   exit;
fi

if [ -z "$8" ]; then
   echo "Specify bound file";
   exit;
fi

if [ -z "$9" ]; then
   echo "Specify ec to reaction file (Kegg/aux/ec-list.txt)";
   exit;
fi

if [ -z "${10}" ]; then
   echo "Specify exchange reactions file";
   exit;
fi

if [ -z "${11}" ]; then
   echo "Specify rxnname file";
   exit;
fi

if [ -z "${12}" ]; then
   echo "Specify pathway file";
   exit;
fi


#This need to be specified to your working dir
CDIR=$5


Kegg=$6
KDIR=$Kegg/kegg-no-general

EQN=$KDIR/nooxygen.eqn

ODIR_PARAM_LABELED=0; #no param labeling. Match with this param in reco-dir and reco-dir.sh

if [ $ODIR_PARAM_LABELED -ne 0 ]; then
    DDIR=$1/reco/$4
else 
    DDIR=$1/reco-$2-$3/$4
fi
mkdir -p $DDIR

#REMOVE_PARTIAL_EC=no
REMOVE_PARTIAL_EC=yes    # remove partial EC numbers from reaction scores

# Convert result into stoichiometric matrix
echo "Generating stoichiometric matrix..."
echo "python $CDIR/network2matrix.py $DDIR $EQN $Kegg/aux/kegg-compounds $DDIR/network"
python $CDIR/network2matrix.py $DDIR $EQN $Kegg/aux/kegg-compounds $DDIR/network
[ $? -ne 0 ] && echo network2matrix.py failed && exit $?

# Convert result into EC list
echo Generating EC list...
echo "python $CDIR/network2eclist.py $DDIR $DDIR/network.ecs"
python $CDIR/network2eclist.py $DDIR $DDIR/network.ecs
[ $? -ne 0 ] && echo network2eclist.py failed && exit $?

# Convert result into EC graph
echo Generating EC graph...
echo "python $CDIR/network2ecgraph.py $DDIR $9"
python $CDIR/network2ecgraph.py $DDIR $9
[ $? -ne 0 ] && echo network2ecgraph.py failed && exit $?

# Convert result into SBML
echo Generating SBML...
taxon=`grep $4 $7|cut -d " " -f 1`
echo $4: taxonomy id $taxon
echo "python $CDIR/network2sbml.py $DDIR $EQN $Kegg/aux/kegg-compounds $taxon $DDIR/$4_$taxon $4_CoReCo $4 $DDIR/$4.sbml"
python $CDIR/network2sbml.py $DDIR $EQN $Kegg/aux/kegg-compounds $taxon $DDIR/$4_$taxon $4_CoReCo $4 $DDIR/$4.sbml ${11} ${12}
[ $? -ne 0 ] && echo network2sbml.py failed && exit $?

echo Post-processing...
echo "python $CDIR/exchange.py $DDIR/$4.sbml ${10} $8"
python $CDIR/exchange.py $DDIR/$4.sbml ${10} $8
[ $? -ne 0 ] && echo exchange.py failed && exit $?
