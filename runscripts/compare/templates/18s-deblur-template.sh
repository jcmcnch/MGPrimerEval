#!/bin/bash

snakemake --configfile config/config-GP13-samples-only.yaml --snakefile Snakefile-compare2tags-18s.smk \
    --cores 20 --use-conda \
    --config cutoff=0.01 pcid=97 denoiser=deblur iLenDeblurTrunc=0 \
    ASVtable=denoiser-output/191004_ASV_info-deblur/EUKs/all-18S-seqs.with-SILVA132-tax.proportions.tsv \
    ASVseqs=denoiser-output/191004_ASV_info-deblur/EUKs/dna-sequences.fasta \
    iLenR1Trunc=220 \
    iLenR2Trunc=180
