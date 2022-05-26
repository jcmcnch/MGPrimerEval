#!/bin/bash -i
conda activate snakemake-env
snakemake --cores 8 --use-conda --snakefile Snakefile-compute.smk --configfile config/tutorial/config.yaml #--conda-create-envs-only #Uncomment out if you want to test conda's ability to set up environments properly
snakemake --cores 8 --use-conda --snakefile Snakefile-classify.smk --configfile config/tutorial/config.yaml
