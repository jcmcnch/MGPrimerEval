## Preamble & Pipeline Rationale

Accurate design of oligonucleotide primers for small subunit ribosomal RNA (SSU rRNA) polymerase chain reaction (PCR) amplicon sequencing (or indeed any PCR-based analysis) determines how quantitiative the resulting data is. So far, primers have been mainly designed based on comprehensive and highly-curated reference databases such as SILVA. This has provided important insights into theoretical primer performance and corrected many flaws, but may be subject to human error/oversight. This pipeline is an attempt to provide an automated, less biased way of evaluating primer performance based on meta-'omics datasets from the environment of interest. It does so with three main modules\*:

* *Compute workflow*: Calculate, from raw (unassembled) metagenome/transcriptome reads, the proportion of reads that perfectly match your primer(s) of interest.
* *Classify workflow*: Using output from the *Compute workflow*, determine which taxa are matched, and which taxa are mismatched. Provide information on primer variants for specific groupings that can allow an investigator to correct biases as appropriate.
* *Compare workflow*: Using output from the *Compute workflow*, compare relative abundances of amplicon sequence variants with those from the *same region* of the 16S/18S molecule. This allows an "apples to apples" comparison between amplicon/metagenomic methods, and can show how well the two methods correspond with one another. 

\* NB: Modules are separated for the sake of clarity, not necessity.

We have described the results of this analysis already for oceanic ecosystems (see preprint [here](https://www.biorxiv.org/content/10.1101/2020.11.09.375543v1)), but we should note that *our pipeline is agnostic to primer/dataset* and thus could be used on any meta-'omics dataset. We hope the instructions below are enough to get you started. If there are any bugs, or questions, please open a github issue above and we'll do our best to help.

## Usage

### Input Requirements:

The only thing you need are raw, *unassembled* paired-end meta'omics data. The pipeline will extract SSU rRNA and analyze them automatically from your data. By default, the pipeline is set to limit the number of retrieved SSU rRNA per sample to 1 million (so it doesn't run too slowly on rRNA-rich data such as metatranscriptomes), but you can tweak this if you'd like to get more data back - just change the `readlimit` parameter in your config file.

## Tutorial

**These steps were tested on a remote server running Ubuntu 16.04 in January 2021.**

I'm going to assume you're familiar with [basic bash command line syntax](https://astrobiomike.github.io/unix/unix-intro), have github installed, and you're using something like `screen` or `tmux` to keep a persistent session alive. I'll also assume you've followed the [snakemake install instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) with the small difference that my conda environment for snakemake is `snakemake-env` not `snakemake`.

### Getting setup: Downloading databases

Our pipeline depends on open-source software, and leans heavily in particular on the wonderful [phyloFlash package](https://github.com/HRGV/phyloFlash/blob/master/README.md) which uses a curated version of the SILVA database to retrieve SSU rRNA from meta'omics data (both 16S and 18S). Before running the tutorial, you will need to download their database. Choose a sensible location, install phyloFlash using conda (the following worked for me: `conda config --set channel_priority true && conda create -n pf sortmerna=2.1b phyloflash && conda config --set channel_priority false`), activate the environment (`conda activate pf`), run their download script in the sensible location (e.g. `phyloFlash_makedb.pl --remote`; this step will take up to a few hours as phyloFlash runs quality-control steps on the database but only needs to be run once), and then make a note of where the database is located so you can add it into your configuration file later.

There is also a splitting step for which you will https://osf.io/e65rs/

### Running the workflow

I usually make a new folder for each new dataset I'm analyzing to keep things organized. For the purposes of this tutorial, let's download the repo into a folder called `MGPrimerEval-tutorial` as follows:

`git clone git@github.com:jcmcnch/MGPrimerEval.git MGPrimerEval-tutorial`

1. Now, `cd` into the folder. If you have data to play with, add it into the folder `intermediate/compute-workflow-intermediate/00-fastq/`. If not, clone [this extremely helpful repo](https://github.com/wwood/ena-fast-download) into the tutorial folder (you should now see a `ena-fast-download` subfolder), install [ascp](https://download.asperasoft.com/download/docs/ascp/3.5.2/html/index.html) (which is proprietary software from IBM that downloads things *very fast*) and use the script `./tutorial/download-BGT.sh` to download a small subset of the [BioGEOTRACES metagenomes](https://www.nature.com/articles/sdata2018176). The script will also put the files into the right spot. *Keep in mind, this is still a fair bit of data! If possible, do it on a work server, not your home network unless you have unlimited bandwidth.*
2. The next step is to set up your configuration file for snakemake to read so it knows what samples to analyze. If you're using the BioGEOTRACES data, there is already a configuration file in `config/tutorial/config.yaml`. If not, you can use the template at `config/config-template.yaml`. *Make sure to modify it to point to database locations and suit your needs*. The most important thing is to add a study name and put your samples at the end (both already done for the tutorial data), put the locations of your databases (i.e. edit the lines starting with `phyloFlashDB` and `bbsplitDBpath`), and the appropriate suffixes for your data (i.e. what should be stripped off to get the sample name; this has also been done for the tutorial config). You can optionally add/remove primers as you see fit (to remove existing primers, just comment them out from the config). 
3. Before you run snakemake, you'll have to activate the conda environment i.e. `conda activate snakemake-env`.
4. Now, run the compute step (upon which the other two modules depend) as follows if you downloaded the tutorial data: `snakemake --cores <# of cores> --use-conda --snakefile Snakefile-compute.smk --configfile config/tutorial/config.yaml`. Snakemake will now automagically install all software dependencies and should begin cranking out data. If you want to test whether everything is working before initiating the run, you can append `-np` to the above command. This will do a dry run and print out some useful output about the steps snakemake plans to execute.

Known issues:

* If you are running into issues with DAG generation (read [snakemake documentation](https://snakemake.readthedocs.io/en/stable/) if you want to know what a DAG is) taking a long time, especially if you have a *lot* of samples, you might need to subset your workflow. This may be fixed in newer versions of snakemake, but was an issue for me about a year ago. You can find some examples of how to do so in the `runscripts` folder.
* bbmap/bbsplit steps sometimes can hang under situations of high RAM use. In this case, just CTRL-C and restart your workflow.

A visual demonstration of the compute module:
![Rule graph for initial (Snakefile-compute.smk) steps](https://github.com/jcmcnch/MGPrimerEval/blob/master/images/Snakefile-compute.svg)
