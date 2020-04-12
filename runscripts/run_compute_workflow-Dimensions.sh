#!/bin/bash
source activate snakemake-env

for item in config/Dimensions_SPOT_CAT/subsetted/*; do 

	snakemake --cores 24 --use-conda --snakefile Snakefile-compute.smk --configfile $item

done

conda deactivate
