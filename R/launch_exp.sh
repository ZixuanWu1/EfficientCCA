#!/bin/bash

# Define the values for the variables

n_values="400"
p_values="100 200 300 400 500"

for n in $n_values; do
for p in $p_values; do
sbatch EfficientCCA/R/exp.sh "$n" "$p"
done
done

