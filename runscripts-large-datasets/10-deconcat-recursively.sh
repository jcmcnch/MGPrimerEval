#!/bin/bash
source activate snakemake-env

for item in `ls subsetted-configs/*yaml`; do

	snakemake --cores 10 --ri --force-use-threads --resources mem_mb=500000 --until deconcat_matches_classifications --configfile $item

done

source deactivate
