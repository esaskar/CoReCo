#$ -S /bin/bash
#$ -N VTTCLusterNetworkReconstruction
#$ -j y
#$ -cwd
#$ -R y


# $1  data dir 
# $2  param-accept
# $3  param-reject
# $4  species
# $5  label-outputdir-with-params
# $8  bound file
# $9  ec-list file


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

#This need to be specified to your working dir
CDIR=$5


Kegg=$6
KDIR=$Kegg/kegg-no-general

EQN=$KDIR/nooxygen.eqn

ERRORDIR=$1/reco-logs/
mkdir -p $ERRORDIR

#ODIR_PARAM_LABELED=no;
ODIR_PARAM_LABELED=yes;


#REMOVE_PARTIAL_EC=no
REMOVE_PARTIAL_EC=yes    # remove partial EC numbers from reaction scores

echo "Reconstructing $4 in $1 (accept $2, reject $3)"
echo "python $CDIR/reco-dir.py $8 $Kegg $KDIR $1 $2 $3 $9 $4 $ODIR_PARAM_LABELED $REMOVE_PARTIAL_EC >$ERRORDIR/$4-$2-$3-out 2>$ERRORDIR/$4-$2-$3-err"
python $CDIR/reco-dir.py $8 $Kegg $KDIR $1 $2 $3 $9 $4 $ODIR_PARAM_LABELED $REMOVE_PARTIAL_EC >$ERRORDIR/$4-$2-$3-out 2>$ERRORDIR/$4-$2-$3-err


[ $? -ne 0 ] && echo reco-dir.py failed && exit $?


