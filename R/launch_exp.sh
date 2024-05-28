#!/bin/bash

# Define the values for the variables

n_values="400"
p_values="200 400 600"

for n in $n_values; do
for p in $p_values; do
sbatch exp.sh "$n" "$p"
done
done
done