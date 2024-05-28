#!/bin/bash

# Define the values for the variables

n_values="200 400 600 800"
p_values="200 300 400"

for n in $n_values; do
for p in $p_values; do
sbatch EfficientCCA/R/exp.sh "$n" "$p"
done
done

