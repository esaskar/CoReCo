#!/bin/bash

mkdir -p fungi-results/knockout
#./perform_knockout_experiment fungi-projects/20120724-species/reco/Scer/network.reactions fungi-results/knockout/Scer
#./perform_knockout_experiment fungi-projects/20120724-species/reco/Scer/network.reactions fungi-results/knockout/Scer-TN
#./perform_knockout_experiment fungi-projects/20120724-species/reco/Scer/network.reactions fungi-results/knockout/Scer-const-delall
#./perform_knockout_experiment fungi-projects/20120724-species/reco/Scer/network.reactions fungi-results/knockout/Scer-const-delcomp
./perform_knockout_experiment fungi-projects/20120724-species/reco/Scer/network.reactions fungi-results/knockout/Scer-noconst-delall
#./perform_knockout_experiment fungi-projects/20120724-species/reco/Scer/network.reactions fungi-results/knockout/Scer-noconst-delcomp
