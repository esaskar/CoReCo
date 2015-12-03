#!/bin/bash

if [[ -z "$5" ]]; then
    echo "usage: $0 name media constraints biomass ddir"
    exit 2;
fi

NAME=$1
MEDIA=$2
CONST=$3
BIOMASS=$4
DDIR=$5

#MEDIA=materials/media/minimal_bicarb
#MEDIA=materials/knockout/media/YPD
#BIOMASS=materials/biomass/iMM904-biomass-20110920 
KEGG=kegg-no-general
TDIR=fungi-results/biomass
mkdir -p $TDIR
TARGET=$TDIR/$NAME.biomass
#DDIR=fungi-projects/20120606/reco-0.02
#CONST="-u materials/gibbs/paula_bounds_040612_ep.txt"
#CONST=
echo "#Output of $0" >$TARGET
echo "#Name: $NAME" >>$TARGET
echo "#Media: $MEDIA" >>$TARGET
echo "#Biomass: $BIOMASS" >>$TARGET
echo "#Constraints: $CONST" >>$TARGET
echo "#Dir: $DDIR" >>$TARGET
for d in $DDIR/*
do
    if [[ -d $d ]]; then
	base=${d##*/};
	echo -n "$base " >>$TARGET;
	echo -n "$base ";
	./fba.py -m $d/network.reactions -i $MEDIA -b $BIOMASS -s $KEGG/nooxygen.eqn -u $CONST 2>/dev/null | grep Yield | cut -c 20- >>$TARGET;
    fi
done
echo
