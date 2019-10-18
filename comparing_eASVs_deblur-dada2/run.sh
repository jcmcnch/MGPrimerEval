#!/bin/bash

mkdir db results

makeblastdb -in 191011_ASV_info-DADA2/191011_GP13_all-16S-seqs.dna-sequences-DADA2.fasta -dbtype nucl -out db/DADA2-eASVS-db

blastn -qcov_hsp_perc 100 -perc_identity 100 -outfmt 6 -query 191004_ASV_info-deblur/191004_GP13_all-16S-seqs.dna-sequences.fasta -db db/DADA2-eASVS-db > results/191018_deblur_eASVs_vs_DADA2.tsv

mkdir results/correlations

./plot-eASVs-against-each-other.py results/191018_deblur_eASVs_vs_DADA2.tsv 191011_ASV_info-DADA2/191011_GP13_all-16S-seqs.with-tax.proportions.DADA2.tsv 191004_ASV_info-deblur/191004_GP13_all-16S-seqs.with-tax.proportions.tsv

mkdir results/graphs

source activate snakemake-env

for item in results/correlations/* ; do ./seaborn-plot-correlations.py $item `basename $item .tsv` results/graphs/`basename $item .tsv`.svg; done
