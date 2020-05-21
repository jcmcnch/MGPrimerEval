#!/bin/bash

snakemake --configfile config/compare/config-GA03-GP13-samples-only.yaml --snakefile Snakefile-compare2tags-16s.smk \
    --cores 20 --use-conda \
    --config cutoff=0.01 pcid=97 denoiser=dada2-new-data-full-table iLenDeblurTrunc=0 \
    ASVtable=config/compare/200514_ASV_info/DADA2/PROKs/all-16S-seqs.with-tax.proportions.tsv \
    ASVseqs=config/compare/200514_ASV_info/DADA2/PROKs/all-16S-seqs.with-tax.proportions.fasta \
    iLenR1Trunc=220 \
    iLenR2Trunc=180 --until plot_ASV_vs_BLAST_results_log_scale #\
#    datestamp='"2019-11-18"' \
#    -np
