#Gene clustering

#! bin/bash

#Run prodigal on trimmed contigs not clustered into vOTUs
input_folder=~/PIMMSgit/data/phages/fasta_files/all_phages/
output_folder=~/PIMMSgit/data/results/gene_clusters/

prodigal -i "$input_folder/all_trimmed.fna" -a "$output_folder/all_genes.faa" -d "$output_folder/all_genes.fna" -p meta -f gff -o "$output_folder/all_genes.gff"