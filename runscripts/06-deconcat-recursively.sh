#!/bin/bash
source activate snakemake-env

snakemake --cores 40 --ri --force-use-threads --resources mem_mb=500000 --until deconcat_classifications

source deactivate
