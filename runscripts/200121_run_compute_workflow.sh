#!/bin/bash
source activate snakemake-env

for item in config/Malaspina/subsetted/*; do 

	echo snakemake --cores 20 --use-conda --snakefile Snakefile-compute.smk --configfile $item

done

conda deactivate
