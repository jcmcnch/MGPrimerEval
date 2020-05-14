#!/bin/bash

snakemake --configfile config/config-GP13-samples-only.yaml --snakefile Snakefile-compare2tags-16s.smk \
    --cores 20 --use-conda \
    --config cutoff=0.01 pcid=100 denoiser=deblur iLenDeblurTrunc=10 \
    ASVtable=denoiser-output/191004_ASV_info-deblur/PROKs/191004_GP13_all-16S-seqs.with-tax.proportions.tsv \
    ASVseqs=denoiser-output/191004_ASV_info-deblur/PROKs/191004_GP13_all-16S-seqs.dna-sequences.fasta 
