#!/bin/bash

snakemake --configfile config/compare/config-GA03-GP13-samples-only.yaml --snakefile Snakefile-compare2tags-16s.smk \
    --cores 20 --use-conda \
    --config cutoff=0.01 pcid=97 denoiser=dada2-old-data iLenDeblurTrunc=0 \
    ASVtable=config/compare/200514_ASV_info/DADA2/PROKs/200519_GA03-GP13_all-16S-seqs.with-tax.proportions.tsv \
    ASVseqs=config/compare/200514_ASV_info/DADA2/PROKs/200519_GA03-GP13_all-16S-seqs.with-tax.proportions.fasta \
    iLenR1Trunc=220 \
    iLenR2Trunc=180 \
    datestamp='"2020-05-20"' \
    --until plot_ASV_vs_BLAST_results_log_scale
