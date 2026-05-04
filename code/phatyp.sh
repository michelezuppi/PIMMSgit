#!/bin/bash

cd ~/miniconda3/envs/PhaTYP/

input_folder=~/PIMMSgit/data/phages/fasta_files
intermediate_folder=~/PIMMSgit/data/results/PHAtyp_results/intermediate
output_folder=~/PIMMSgit/data/results
code_folder=~/miniconda3/envs/phatyp/bin

{ time \
  taskset -c 16-30 python "$code_folder/preprocessing.py" \
    --contigs "$input_folder/all_species_representative.fna" \
    --len 1500 \
    --midfolder "$intermediate_folder" && \
  taskset -c 0-15 python "$code_folder/PhaTYP.py" \
    --out "$output_folder/PHAtyp_results/phatyp_output.csv" \
    --midfolder "$intermediate_folder" ; \
} 2> "$output_folder/PHAtyp_results/phatyp_run_time.txt"

