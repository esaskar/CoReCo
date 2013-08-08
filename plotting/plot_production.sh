#!/bin/bash

cut -f 1,2,4 materials/growth/20111212b |grep -v NA >materials/growth/20111212b.frac.cut
cut -f 1,2,3 materials/growth/20111212b |grep -v NA >materials/growth/20111212b.abs.cut
./to_matrix.py materials/growth/20111212b.frac.cut materials/growth/20111212b.frac.matrix
./to_matrix.py materials/growth/20111212b.abs.cut materials/growth/20111212b.abs.matrix
#grep -v NA materials/growth/20111212b.matrix >materials/growth/20111212b.matrix2
echo "Writing to materials/growth/20111212b.pdf"
R --no-save <plot_production.R
