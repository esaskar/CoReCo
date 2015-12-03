#!/bin/bash

MEDIA=materials/media/minimal_bicarb
BIOMASS=materials/chemical-production/top_value_added_chemicals
KEGG=kegg-no-general
TDIR=fungi-results/chemical/noconst
mkdir -p $TDIR
TARGET=$TDIR
DDIR=fungi-projects/20120724-species/reco
#CONST="-u materials/gibbs/paula_bounds_040612_ep.txt"
CONST=

# Calculate yields for whole KEGG
./fba.py -m kegg-no-general/rlist -s $KEGG/nooxygen.eqn -i $MEDIA -b $BIOMASS $CONST -c >$TARGET/KEGG

for d in $DDIR/*
do
    if [[ -d $d ]]; then
	base=${d##*/};
	echo -n "$base ";
	echo "#Output of $0" >$TARGET/$base
	echo "#Species: $base" >>$TARGET/$base
	echo "#Media: $MEDIA" >>$TARGET/$base
	echo "#Biomass: $BIOMASS" >>$TARGET/$base
	echo "#Constraints: $CONST" >>$TARGET/$base
	echo "#Dir: $DDIR" >>$TARGET/$base
	./fba.py -m $d/network.reactions -i $MEDIA -b $BIOMASS -s $KEGG/nooxygen.eqn $CONST -c 2>/dev/null >>$TARGET/$base;
#	./fba.py -m $d/network.reactions -i $MEDIA -b $BIOMASS -s $KEGG/nooxygen.eqn $CONST -c >>$TARGET/$base;
    fi
done
echo

./fungi_plot_biomass_single.sh $TARGET fungi-results/chemical/noconst.fluxes
