#!/bin/bash
source activate snakemake-env

snakemake --use-conda --cores 40 --ri --force-use-threads --resources mem_mb=500000 --until classify_matches_subsample --configfile config.yaml

source deactivate
