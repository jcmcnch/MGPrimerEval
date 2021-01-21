## Preamble & Pipeline Rationale

Accurate design of oligonucleotide primers for small subunit ribosomal RNA (SSU rRNA) polymerase chain reaction (PCR) amplicon sequencing (or indeed any PCR-based analysis) determines how quantitiative the resulting data is. So far, primers have been mainly designed based on comprehensive and highly-curated reference databases such as SILVA. This has provided important insights into theoretical primer performance and corrected many flaws, but may be subject to human error/oversight. This pipeline is an attempt to provide an automated, less biased way of evaluating primer performance based on meta-'omics datasets from the environment of interest. It does so with three main modules\*:

* *Compute workflow*: Calculate, from raw (unassembled) metagenome/transcriptome reads, the proportion of reads that perfectly match your primer(s) of interest.
* *Classify workflow*: Using output from the *Compute workflow*, determine which taxa are matched, and which taxa are mismatched. Provide information on primer variants for specific groupings that can allow an investigator to correct biases as appropriate.
* *Compare workflow*: Using output from the *Compute workflow*, compare relative abundances of amplicon sequence variants with those from the *same region* of the 16S/18S molecule. This allows an "apples to apples" comparison between amplicon/metagenomic methods, and can show how well the two methods correspond with one another. 

\* NB: Modules are separated for the sake of clarity, not necessity.

We have described the results of this analysis already for oceanic ecosystems (see preprint [here](https://www.biorxiv.org/content/10.1101/2020.11.09.375543v1)), but we should note that *our pipeline is agnostic to primer/dataset* and thus could be used on any meta-'omics dataset. We hope the instructions below are enough to get you started. If there are any bugs, or questions, please open a github issue above and we'll do our best to help.

## Usage

### Input Requirements:

The only thing you need are raw, *unassembled* paired-end meta'omics data. The compute step by default limits the number of retrieved SSU rRNA per sample to 1 million, but you can tweak this if you'd like to get more data back.

### Tutorial

I'm going to assume you're familiar with [basic bash command line syntax](https://astrobiomike.github.io/unix/unix-intro), have github installed, and you're using something like screen or tmux to keep a persistent session alive. I'll also assume you've followed the [snakemake install instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) with the small difference that my conda environment for snakemake is `snakemake-env` not `snakemake` as found in the install instructions.

I usually make a new folder for each new dataset I'm analyzing to keep things organized. For the purposes of this tutorial, let's download the repo into a folder called `MGPrimerEval-tutorial` as follows:

`git clone git@github.com:jcmcnch/MGPrimerEval.git MGPrimerEval-tutorial`

1. Now, once you `cd` into the folder, everything is set up for you. You just need to add data into the folder `intermediate/compute-workflow-intermediate/00-fastq/`. You can use your own data, or if you'd like you can clone [this extremely helpful repo](https://github.com/wwood/ena-fast-download) into the tutorial folder, install [ascp](https://download.asperasoft.com/download/docs/ascp/3.5.2/html/index.html) and use the template script `tutorial/download-BGT.sh` to download a small subset of the [BioGEOTRACES metagenomes](https://www.nature.com/articles/sdata2018176). *Keep in mind, this is still a lot of data! Maybe do it on a work server, not your home network unless you have unlimited bandwidth.*
2. Create a configuration file from a template. One is provided for this tutorial, but you might do it as follows: bash munging. *Make sure to modify it to suit your needs*. You can add/remove primers as you see fit, and you should put in the appropriate suffixes.
3. Install Snakemake (give link; note has not been tested on Windows)
4. Activate environment (assume snakemake-env)
5. Now you can run it simply with the following commands:
6. If you are running into issues with DAG generation taking a long time.

Known issues:

* bbmap/bbsplit steps sometimes can hang under situations of high RAM use. In this case, just CTRL-C and restart your workflow.

This is a snakemake pipeline that determines primer mismatches to SSU rRNA in metagenomic data. It can be customized with whatever primer sequence you want to put in (just change the config.yaml file), you just need to make sure to reverse complement the reverse primer. You'll probably also have to customize the "00-fastq" file input unless you have the exact same input format.

![Rule graph for initial (Snakefile-compute.smk) steps](https://github.com/jcmcnch/MGPrimerEval/blob/master/images/Snakefile-compute.svg)
