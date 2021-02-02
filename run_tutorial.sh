#!/bin/bash -i
conda activate snakemake-env
snakemake --cores 80 --use-conda --snakefile Snakefile-compute.smk --configfile config/tutorial/config.yaml
