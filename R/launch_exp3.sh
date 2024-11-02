#!/bin/bash

# Define the values for the variables

n_values="400"
p_values="200 400 600 800 1000"
s_values="5 10 15 20 25"

for n in $n_values; do
for p in $p_values; do
for s in $s_values; do
sbatch EfficientCCA/R/exp3.sh "$n" "$p" "$s"
done
done
done

