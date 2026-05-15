#Run MMseq2


#!bin/bash

input_folder=~/PIMMSgit/data/results/gene_clusters
mkdir "$input_folder/cluster"
mkdir "$input_folder/msa_db"
mkdir "$input_folder/profile_db"
mkdir "$input_folder/consensus_db"
mkdir "$input_folder/profile_consensus_search_db"
mkdir "$input_folder/profile_cluster_db"
mkdir "$input_folder/merged_cluster_db/"
mkdir "$input_folder/nr_output"


#Database creation

#mmseqs createdb "$input_folder/all_genes.faa" "$input_folder/all_genes_DB"

## Step 1 — Full clustering
mmseqs cluster "$input_folder/all_genes_DB" "$input_folder/cluster/all_genes_DB_clu" "$input_folder/tmp" \
    -s 7.0 -c 0.8 --cov-mode 0 --cluster-steps 3 --cluster-reassign 1 --kmer-per-seq 50 \
    --threads 16
    
#Reset the tmp folder   
rm -rf "$input_folder/tmp"

# Generate MSA for each cluster using centre star alignment
mmseqs result2msa "$input_folder/all_genes_DB" "$input_folder/all_genes_DB" \
    "$input_folder/cluster/all_genes_DB_clu" "$input_folder/msa_db/msa_db" \
    --msa-format-mode 2 --threads 16

# Get consensus sequences from profiles
mmseqs msa2profile --msa-type 2 --match-mode 1 --match-ratio 0.5 "$input_folder/msa_db/msa_db" "$input_folder/profile_db/profile_db" --threads 16

#Get consensus sequence
mmseqs profile2consensus "$input_folder/profile_db/profile_db" "$input_folder/consensus_db/consensus_db" --threads 16

# Profile-vs-protein search

mmseqs search --cov-mode 0 -c 0.9 -s 8.0 -e 1e-4 --add-self-matches 1 -a 1 "$input_folder/profile_db/profile_db" \
  "$input_folder/consensus_db/consensus_db" "$input_folder/profile_consensus_search_db/profile_consensus_search_db" "$input_folder/tmp" --threads 16

#Reset the tmp folder   
rm -rf "$input_folder/tmp"

# Step 2 — Cluster original proteins guided by search results
mmseqs clust "$input_folder/profile_db/profile_db" "$input_folder/profile_consensus_search_db/profile_consensus_search_db" \
  "$input_folder/profile_cluster_db/profile_cluster_db" --threads 16

# Merge the clusters generated

mmseqs mergeclusters "$input_folder/all_genes_DB" "$input_folder/merged_cluster_db/merged_cluster_db" "$input_folder/cluster/all_genes_DB_clu" \
  "$input_folder/profile_cluster_db/profile_cluster_db" --threads 16

# Final outputs
mmseqs createtsv "$input_folder/all_genes_DB" "$input_folder/all_genes_DB" "$input_folder/merged_cluster_db/merged_cluster_db"  \
  "$input_folder/nr_output/protein_clusters.tsv"