#!/bin/bash
source activate snakemake-env

for item in config/BioGEOTRACES/subsetted-illinois/*yaml; do 

	snakemake --cores 56 --use-conda --snakefile Snakefile-compute.smk --configfile $item

done

conda deactivate
