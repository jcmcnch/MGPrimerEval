#!/bin/bash
source activate snakemake-env

for item in `ls subsetted-configs/*yaml`; do

	snakemake --cores 40 --ri --force-use-threads --resources mem_mb=500000 --until deconcat_classifications --configfile $item

done

source deactivate