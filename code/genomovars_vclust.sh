#!/bin/bash
input_folder=~/PIMMSgit/data/phages/unified/vclust
output_folder=~/PIMMSgit/data/phages/unified/vclust

{ time \
  taskset -c 0-15 vclust cluster \
    -i "$input_folder/_ani.tsv" \
    -o "$output_folder/_identical.tsv" \
    --ids "$input_folder/_ani.ids.tsv" \
    --algorithm cd-hit \
    --metric ani \
    --ani 0.995 \
    --qcov 0.85 ; \
} 2> "$output_folder/vclust_run_time.txt"


