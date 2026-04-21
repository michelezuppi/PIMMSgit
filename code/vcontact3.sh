#!/bin/bash
#paths for input and output folders
input_folder=~/PIMMSgit/data/
output_folder=~/PIMMSgit/data/results/
#number of threads
threads=16
#runtime
start=$SECONDS

vcontact3 run \
  --nucleotide "$input_folder/phages/all_species_representatives.fna" \
  --output "$output_folder/vcontact3" \
  --db-domain "prokaryotes" \
  --db-version 230 \
  --db-path "$input_folder/databases/virus-host-db" \
  --threads "$threads"

elapsed=$((SECONDS - start ))
echo "Runtime: $(( elapsed / 60 )) minutes and $(( elapsed % 60 )) seconds" > "$output_folder/vcontact3/runtime.txt"
