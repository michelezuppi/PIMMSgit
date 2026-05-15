#!/bin/bash
input_folder=~/PIMMSgit/data/results/gene_clusters
mkdir -p "$input_folder/sequence_identity"

# Get representative sequences as MMseqs2 DB
 mmseqs createdb "$input_folder/merged_cluster_db/merged_cluster_db" \
     "$input_folder/all_genes_DB" 

#Create subdatabase
    mmseqs createsubdb \
        "$input_folder/merged_cluster_db/merged_cluster_db" \
        "$input_folder/all_genes_DB" \
        "$input_folder/merged_cluster_db/merged_cluster_db_rep"

# Convert to FASTA
 mmseqs convert2fasta \
     "$input_folder/merged_cluster_db/merged_cluster_db_rep" \
     "$input_folder/nr_output/merged_cluster_rep.faa"

# Step 1 — Create fresh standalone DB from FASTA (for identity clustering)
mmseqs createdb \
    "$input_folder/nr_output/merged_cluster_rep.faa" \
    "$input_folder/sequence_identity/merged_cluster_db_rep"

#Do the clustering

for identity in 0.90 0.70 0.50 0.30; do

    outdir="$input_folder/sequence_identity/identity_${identity}"
    mkdir -p "$outdir"
    mkdir -p "$input_folder/tmp_${identity}"

    echo "Running clustering at identity ${identity}..."

    mmseqs cluster \
        "$input_folder/sequence_identity/merged_cluster_db_rep" \
        "$outdir/cluster_${identity}" \
        "$input_folder/tmp_${identity}" \
        --min-seq-id ${identity} \
        --cluster-mode 0 \
        --cluster-steps 3 \
        --cluster-reassign 1 \
        --alignment-mode 3 \
        --threads 16

    mmseqs createtsv \
        "$input_folder/sequence_identity/merged_cluster_db_rep" \
        "$input_folder/sequence_identity/merged_cluster_db_rep" \
        "$outdir/cluster_${identity}" \
        "$outdir/cluster_${identity}.tsv"

    echo "Done with identity ${identity}. Clusters: $(cut -f1 $outdir/cluster_${identity}.tsv | sort -u | wc -l)"

    rm -rf "$input_folder/tmp_${identity}"

done 2> "$input_folder/READerror.txt"