#!/bin/bash
#
# This script performs biomass experiments for the fungi reconstruction paper
# 
# 18.7.2012 / EP


MEDIA=/users/jianhou/work/model/model_training/reconstruction/plotting/minimal_bicarb
CONST=/users/jianhou/work/model/model_training/reconstruction/plotting/paula_bounds_040612_ep.txt
CONST_CO2=materials/gibbs/paula_bounds_040612_ep_CO2assim.txt
CONST_NA=-
BIOMASS=/users/jianhou/work/model/model_training/reconstruction/plotting/iMM904-biomass-20110920 
DDIR=/users/jianhou/work/model/model_training/reconstruction/project/reco/Scer

TDIR=/users/jianhou/work/model/model_training/reconstruction/project/reco

NAME=minimal_bicarb-iMM904bm-0.02-const
echo "Biomass / $NAME..."
#./fungi_check_biomass.sh $NAME $MEDIA $CONST $BIOMASS $DDIR
echo "Biomass components / $NAME..."
#./fungi_check_biomass_comps.sh $NAME $MEDIA $CONST $BIOMASS $DDIR
./fungi_plot_biomass_single.sh $TDIR/$NAME $TDIR/$NAME.bmcomp.fluxes

NAME=minimal_bicarb-iMM904bm-0.02-co2const
echo "Biomass / $NAME..."
#./fungi_check_biomass.sh $NAME $MEDIA $CONST_CO2 $BIOMASS $DDIR
echo "Biomass components / $NAME..."
#./fungi_check_biomass_comps.sh $NAME $MEDIA $CONST_CO2 $BIOMASS $DDIR
./fungi_plot_biomass_single.sh $TDIR/$NAME $TDIR/$NAME.bmcomp.fluxes

NAME=minimal_bicarb-iMM904bm-0.02-noconst
echo "Biomass / $NAME..."
#./fungi_check_biomass.sh $NAME $MEDIA $CONST_NA $BIOMASS $DDIR
echo "Biomass components / $NAME..."
#./fungi_check_biomass_comps.sh $NAME $MEDIA $CONST_NA $BIOMASS $DDIR
./fungi_plot_biomass_single.sh $TDIR/$NAME $TDIR/$NAME.bmcomp.fluxes

