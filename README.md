This is a snakemake pipeline that determines primer mismatches to SSU rRNA in metagenomic data. It can be customized with whatever primer sequence you want to put in (just change the config.yaml file), you just need to make sure to reverse complement the reverse primer. You'll probably also have to customize the "00-fastq" file input unless you have the exact same input format.

![Rule graph for initial (Snakefile-compute.smk) steps](https://github.com/jcmcnch/MGPrimerEval/blob/master/images/Snakefile-compute.svg)
