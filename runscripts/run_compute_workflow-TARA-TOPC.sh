#!/bin/bash -i
conda activate snakemake-env

for item in config/TARA-TOPC/subsetted/*; do 

	snakemake --cores 60 --use-conda --snakefile Snakefile-compute.smk --configfile $item

done

conda deactivate
