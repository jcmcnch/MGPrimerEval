#!/bin/bash

snakemake --configfile config/config-GP13-samples-only.yaml --snakefile Snakefile-compare2tags-18s.smk \
    --cores 20 --use-conda --latency-wait 120 \
    --config cutoff=0.01 pcid=100 denoiser=dada2 \
    ASVtable=denoiser-output/191011_ASV_info-DADA2/EUKs/all-18S-seqs.with-SILVA132-tax.proportions.tsv \
    ASVseqs=denoiser-output/191011_ASV_info-DADA2/EUKs/all-18S-seqs.with-SILVA132-tax.dna-sequences.fasta \
    iLenR1Trunc=220 \
    iLenR2Trunc=180
