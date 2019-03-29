#!/bin/bash
source activate snakemake-env

for item in `ls subsetted-configs/*yaml`; do

	echo snakemake --cores 80 --ri --force-use-threads --resources mem_mb=500000 --until compute_percentages --configfile $item 

done

conda deactivate
