#!/bin/bash

# Define the values for the variables

n_values="400"
p_values="200 300 400 500 600"
s_values="5 10 15"

for n in $n_values; do
for p in $p_values; do
for s in $s_values; do
sbatch EfficientCCA/R/exp.sh "$n" "$p" "$s"
done
done
done

