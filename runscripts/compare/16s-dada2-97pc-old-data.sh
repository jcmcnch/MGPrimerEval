#!/bin/bash

snakemake --configfile config/compare/config-GP13-samples-only.yaml --snakefile Snakefile-compare2tags-16s.smk \
    --cores 20 --use-conda \
    --config cutoff=0.01 pcid=97 denoiser=dada2-old-data iLenDeblurTrunc=0 \
    ASVtable=denoiser-output/191011_ASV_info-DADA2/PROKs/191011_GP13_all-16S-seqs.with-tax.proportions.DADA2.tsv \
    ASVseqs=denoiser-output/191011_ASV_info-DADA2/PROKs/191011_GP13_all-16S-seqs.dna-sequences-DADA2.fasta \
    iLenR1Trunc=220 \
    iLenR2Trunc=180 --until plot_ASV_vs_BLAST_results_log_scale #\
#    datestamp='"2019-11-18"' \
#    -np
