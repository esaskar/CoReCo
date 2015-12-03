#!/bin/bash

cut -f 1,2,4 materials/chemical-production/20111212.yields.filtered >materials/chemical-production/20111212.yields.cut
./to_matrix.py materials/chemical-production/20111212.yields.cut materials/chemical-production/20111212.yields.matrix
grep -v NA materials/chemical-production/20111212.yields.matrix >materials/chemical-production/20111212.yields.matrix2
echo "Writing to materials/chemical-production/20111212.yields.png"
R --no-save <plot_chemical_prod.R
