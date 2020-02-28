#!/bin/bash
source activate snakemake-env

for item in config/TARA-gt3uM/subsetted/*; do 

	snakemake --cores 24 --use-conda --snakefile Snakefile-compute.smk --configfile $item

done

conda deactivate
