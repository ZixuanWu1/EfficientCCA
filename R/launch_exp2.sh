#!/bin/bash

# Define the values for the variables

n_values="400"
p_values="200"
s_values="5 10 15"
r_values="2 3 4 5"

for n in $n_values; do
for p in $p_values; do
for s in $s_values; do
for r in $r_values; do
sbatch EfficientCCA/R/exp2.sh "$n" "$p" "$s" "$r"
done
done
done
done

