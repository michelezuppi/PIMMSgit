#!/bin/bash
#paths for input and output folders
input_folder=~/PIMMSgit/data/
output_folder=~/PIMMSgit/data/results/
#number of threads
threads=16

{ time vcontact3 run \
  --nucleotide "$input_folder/all_species_representatives.fna" \
  --output "$output_folder/vcontact3" \
  --db-domain "prokaryotes" \
  --db-version 228 \
  --db-path "$input_folder/databases/virus-host-db" \
  --threads "$threads" ; } 2> "$output_folder/vcontact3/runtime.txt"
