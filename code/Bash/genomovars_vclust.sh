#!/bin/bash
input_folder=~/PIMMSgit/data/phages/unified/vclust
output_folder=~/PIMMSgit/data/phages/unified/vclust

{ time \
  taskset -c 0-15 vclust cluster \
    -i "$input_folder/_ani.tsv" \
    -o "$output_folder/_genomovars.tsv" \
    --ids "$input_folder/_ani.ids.tsv" \
    --algorithm complete \
    --metric tani \
    --tani 0.995 \
    --out-repr ; \
} 2> "$output_folder/vclust_run_time.txt"


