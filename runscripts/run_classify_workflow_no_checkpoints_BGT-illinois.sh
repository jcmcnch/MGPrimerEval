#!/bin/bash
source activate snakemake-env

cores=24

#Run first part of the workflow but use subsetted-illinois configs to speed up DAG generation and execution (particularly important for large datasets)

for item in config/BioGEOTRACES/subsetted-illinois/*; do

        snakemake --snakefile Snakefile-classify.smk  --configfile $item --use-conda --cores $cores --until classify_mismatches

done

for item in config/BioGEOTRACES/subsetted-illinois/*; do

        snakemake --snakefile Snakefile-classify.smk  --configfile $item --use-conda --cores $cores --until classify_matches_subsample

done

for item in config/BioGEOTRACES/subsetted-illinois/*; do

	snakemake --snakefile Snakefile-classify.smk  --configfile $item --use-conda --cores $cores --until classify_cyano_fraction_phytoRef

done

#Then can run the subsequent steps with a master config file or subsetted-illinois as before - need to decide

snakemake --snakefile Snakefile-classify.smk  --configfile config/BioGEOTRACES/config-BioGEOTRACES-master-illinois.yaml --use-conda --cores $cores --until make_tax_matchVSmismatch_barplots

conda deactivate
