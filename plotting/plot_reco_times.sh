#!/bin/bash

sort reco_times >reco_times.sorted
grep -E -v "N[[:digit:]]+" reco_times.sorted >reco_times.species
R --no-save <plot_reco_time.R
