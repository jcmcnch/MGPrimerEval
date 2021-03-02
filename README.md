Just want to run the pipeline? [Jump to usage instructions.](https://github.com/jcmcnch/MGPrimerEval#running-the-pipeline-with-your-own-data)

## Preamble
Accurate design of oligonucleotide primers for small subunit ribosomal RNA (SSU rRNA) polymerase chain reaction (PCR) amplicon sequencing (or indeed any PCR-based analysis) determines how quantitiative the resulting data is. So far, primers have been mainly designed based on comprehensive and highly-curated reference databases such as SILVA. This has provided important insights into theoretical primer performance and corrected many flaws, but may be subject to human error/oversight. This pipeline is an attempt to provide an automated, less biased way of evaluating primer performance based on meta-'omics datasets from the environment of interest. It does so with three main modules\*:

* *Compute workflow*: Calculate, from raw (unassembled) metagenome/transcriptome reads, the proportion of reads that perfectly match your primer(s) of interest.
* *Classify workflow*: Using output from the *Compute workflow*, determine which taxa are matched, and which taxa are mismatched. Provide information on primer variants for specific groupings that can allow an investigator to correct biases as appropriate.
* *Compare workflow*: Using output from the *Compute workflow*, compare relative abundances of amplicon sequence variants with those from the *same region* of the 16S/18S molecule. This allows an "apples to apples" comparison between amplicon/metagenomic methods, and can show how well the two methods correspond with one another. 

\* NB: Modules are separated for the sake of clarity, not necessity.

We have described the results of this analysis already for oceanic ecosystems (see preprint [here](https://www.biorxiv.org/content/10.1101/2020.11.09.375543v1)), but we should note that *our pipeline is agnostic to primer/dataset* and thus could be used on any meta-'omics dataset. We hope the instructions below are enough to get you started. If there are any bugs, or questions, please open a github issue above and we'll do our best to help.

## Detailed Overview

### Motivating Scientific Questions

This pipeline is designed to address several related questions at different levels of detail:
1. What fraction of environmental SSU rRNA fragments match oligonucleotide primer sequences\* in a given environment at 0, 1, and 2-mismatch thresholds? *i.e., how well do primers theoretically perform for a given environment / dataset?*
2. What taxa are perfectly matched by my primers? What taxa are not, and therefore likely to be inaccurately quantified using PCR amplicon barcoding methods? *i.e. what are the taxonomic blindspots of a given oligo / primer set?*
3. How can I improve a given oligonucleotide primer to improve its performance on a given dataset/environment? *i.e. can I create an "optimal" primer set for my environment?*
4. How well do results from PCR primers compare with metagenomic reads? *i.e. are the two methods consistent or inconsistent with one another?*

\*Note: the pipeline *does* handle degenerate primers as a query, and returns a perfect match if one of the variants specified in your degenerate sequence matches the target.

### Input Requirements:

The only thing you need are raw, *unassembled* paired-end meta'omics data such as metagenomes or metatranscriptomes. They should be compressed in gzip format (suffix=gz). Merged read pairs or single-end reads are not currently supported. It is important that data have not been filtered or assembled, since the goal of this pipeline is to recover the underlying environmental patterns with respect to primer matches/mismatches (which assumes that your metagenome/-transcriptome is an accurate representation of the environment in question). 

### Recommended operating systems

The recommended operating system to run the pipeline is Linux. The pipeline has been tested on Ubuntu 16.04, Ubuntu 20.04, and CentOS. It may work on Mac/Windows so long as you can install snakemake and conda, but has not been tested on these systems.

### Overview of Pipeline Steps:

The pipeline steps are roughly as follows:
1. Extract SSU rRNA fragments (16S and 18S) from your input data\*
2. Quality control retrieved fragments
3. Split SSU rRNA fragments into 4 categories:
- *Archaea*
- *Eukarya* (18S)
- *Bacteria* (not including *Cyanobacteria*)
- *Cyanobacteria* (including chloroplast 16S sequences derived from eukaryotic phytoplankton)
4. Align to a reference and subset these fragments to SSU rRNA region
5. Compare PCR primers specified in the config file with your SSU rRNA fragments
6. Classify matching/mismatching fragments to get their taxonomy
7. Summarize data

\*By default, the pipeline is set to limit the number of retrieved SSU rRNA per sample to 1 million (so it doesn't run too slowly on rRNA-rich data such as metatranscriptomes), but you can tweak this if you'd like to get more data back - just change the `readlimit` parameter in your config file.

### Expected output files

*NB: By default, the pipeline keeps all intermediate files except for the fastp processed raw reads, but you can change this behaviour by putting `temp()` around any output files you wish to discard. That being said, the processed data files should be considerably smaller than your raw data. I also have a cleanup script in the repository you can use to compress and remove some unnecessary intermediates if you're running out of space (`scripts/compress-cleanup-MGPrimerEval.sh`).*

If you run just the *compute* workflow, you will get:
- phyloFlash summaries of taxa present in your metagenome/-transcriptome (can be useful to make sure your sample is what you think it is)
- QC'd SSU rRNA fragments
- Sorted SSU rRNA fragments (into categories noted above)
- Aligned SSU rRNA fragments
- Fragments subsetted from alignments to primer regions, and sorted into matching and non-matching at 0, 1, and 2-mismatch thresholds

If you additionally run the *classify* workflow using output from above, you will get:
- A graphical summary of the overall results
- Taxonomic assignments for above matching / mismatching fragments (at the order level)
- Graphical summaries of the taxonomic matches / mismatches
- Tabular summaries across your whole dataset, indicating the proportions of order-level taxonomic groups are matched or mismatched for each primer/group/mismatch threshold (e.g. summarizing mismatches to *Archaea* for the primer 515Y at a 0-mismatch threshold)
- Tabular summaries of the primer variants, both for the groups indicated above *and* for taxa that have abundant mismatches (so you can correct these mismatches or at least know what they are, even if the taxa are rare)

The *compare* workflow (only if you have paired metagenomes and amplicon sequences you want to intercompare):
- SSU rRNA fragments subsetted to the primer region (done separately for 16S and 18S)
- A BLASTn-based comparison between MG SSU rRNA fragments and amplicon sequence variants (using ASVs as a BLAST database and the MG SSU rRNA as query)
- A direct intercomparison between taxonomic groups found in MG SSU rRNA and ASVs *from the same sample*, summarized in graphical and tabular format (includes R^2 values of relative abundances; see manuscript text for more details)

### Open-source software and database dependencies

This pipeline depends on a number of amazing free and open source software packages such as:

- [Snakemake](https://snakemake.readthedocs.io/en/stable/) which runs the whole workflow.
- The wonderful [phyloFlash package](https://github.com/HRGV/phyloFlash/blob/master/README.md) which uses a curated version of the SILVA database to retrieve SSU rRNA from meta'omics data (both 16S and 18S).
- [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/), the "Swiss Army Knife" of bioinformatics - which we use for splitting the SSU rRNA into different categories as well as data manipulation (bbtools is also a dependency of phyloFlash).
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/) which is used for matching the primers to metagenomic reads.
- [VSEARCH](https://github.com/torognes/vsearch) which is used to taxonomically classify reads.
- [fastp](https://github.com/OpenGene/fastp) which is used to quality-control raw reads.
- [seqtk](https://github.com/lh3/seqtk) which is used for basic sequence manipulation.
- [pyNAST](https://github.com/biocore/pynast) which is used to align SSU rRNA sequences against a reference.

Our taxonomic classification and splitting steps also heavily depend on the [SILVA database project](https://www.arb-silva.de/), which is an expert-curated database that is the most comprehensive SSU rRNA database for diverse global environments.

## Running the pipeline with your own data

The following are instructions to get the pipeline set up for your own datasets. There are also [instructions below for downloading and process example data](https://github.com/jcmcnch/MGPrimerEval#example-data-eg-if-you-just-want-to-testverify-the-functionality-of-the-pipeline-on-your-system) if you just want to test the mechanics and make sure it runs on your system. You will still need to follow most of the setup below, with the exception of adding your sample names to the config (a config is already provided with the sample names for these example data).

These instructions assume you have familiarity with [basic bash command line syntax](https://astrobiomike.github.io/unix/unix-intro), have github installed, and you're using something like `screen` or `tmux` to keep a persistent session alive for remote servers. 

### Setting up and activating a snakemake conda environment

Assuming you have the [python3 version of miniconda](https://conda.io/en/latest/miniconda.html) installed, install snakemake into its own environment and activate it as follows:

```
#install mamba, which is faster than conda
conda install -c conda-forge mamba

#create an environment named snakemake-env
mamba create -c conda-forge -c bioconda -n snakemake-env snakemake
conda activate snakemake-env
```

[Source for install instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

### Cloning the repository and adding raw data

Now, clone the repo into a new folder we'll call `myDataset` and enter that folder.

```
git clone https://github.com/jcmcnch/MGPrimerEval.git myDataset
cd myDataset
```

Link your raw data into the input folder (the `ln -s` "softlink" prevents data duplication). For example:

`ln -s /full/path/to/your/data/*gz intermediate/compute-workflow-intermediate/00-fastq/`

### Downloading databases for phyloFlash, SSU rRNA splitting, SSU rRNA classification with VSEARC, and adding `uclust` to your path

1. PhyloFlash Database for Retrieving SSU rRNA fragments

Install phyloFlash into its own environment and run the database install script (will take some hours as phyloFlash runs quality-control steps on the database but this only needs to be run once):

```
#Install mamba into your base environment (if not already installed)
conda install -c conda-forge mamba

#Use mamba to create phyloFlash environment
mamba create -c conda-forge -c bioconda --name pf sortmerna=2.1b phyloflash

#If you're getting errors, you may need to run `conda update conda` or do a fresh install of miniconda if updating is not easy (sometimes you get all sorts of incompatibilities which can just be solved by a fresh install)
conda activate pf

#change directory to suit your needs
mkdir -p ~/databases/phyloFlash-db/ ; cd ~/databases/phyloFlash-db/

#run database download/construction script
#will take a few hours while the database is downloaded and QC'd
phyloFlash_makedb.pl --remote
```

\*Note: `sortmerna` is not used in this pipeline but is part of the conda recipe.

2. Database for Splitting SSU rRNA fragments

If you do not already have bbmap installed locally, make a conda environment called bbmap-env:

`mamba create -c bioconda --name bbmap-env bbmap`

Now, download and create the database (should only take a few minutes):

```
#activate environment
conda activate bbmap-env

#make directory, enter it
mkdir -p ~/databases/bbsplit-db/ ; cd ~/databases/bbsplit-db

#download database files from OSF
for item in kv3xp eux4r npb2k 4qtev s5j6q 5jmkv eahds ; do curl -O -J -L https://osf.io/$item/download ; done

#make databases
chmod u+x make-dbs-bbsplit.sh ; ./make-dbs-bbsplit.sh
```

3. Databases for classifying SSU rRNA fragments with VSEARCH

**NB: If these databases are not set up correctly, the *Classify* workflow will run but will generate empty output files.** 

The following code will download the `udb`-formatted files you can use for the *Classify* workflow (not needed if you only want to run the *Compute* workflow):

```
#make directory, enter it
mkdir -p ~/databases/VSEARCH_db/ ; cd ~/databases/VSEARCH_db

#download database files from OSF
for item in  IDS ; do curl -O -J -L https://osf.io/$item/download ; done

#get full paths for config file
ls $PWD/*udb
```

4. Getting the `uclust` executable:

**Please note that the pipeline will still run if `uclust` is not set up correctly, but will produce empty output files for the alignment step (meaning you won't get any results). So double-check that the path you provide in the config below is accurate.**

The alignment steps in this pipeline currently depend on `pyNAST`, which also depends on `uclust`. However, the `uclust` executable is not available through standard repositories as it is not open-source. You may have access to `uclust` (e.g. from an older install of qiime), but you can also just email me at mcnichol at alum dot mit dot edu and I'll send you the binary. I have [been given permission](https://github.com/biocore/pynast/issues/21) to distribute the executable I used by email by the author of `uclust`.

All you need to do is put the binary in a sensible location, and make a note of the full path to add to the config file (next section). This way, when you run the workflow, this executable will be found. For example:

```
#enter your analysis folder
cd /home/jesse/MGPrimerEval-tutorial

#make a directory for binaries, enter directory
mkdir bin ; cd bin

#copy uclust into this directory
cp /location/of/uclust .

#make executable
chmod a+x uclust

#get full path to put into the config file (below)
pwd
```

### Setting up your configuration file

The template configuration file comes pre-set with a number of primers that we tested in our study. If you just want to test these primers on your samples, all you have to do is add your sample names at the end in the format `  sample : sample`. I suggest making a new folder and config file for your analysis to keep things organized. If your forward reads end with `_1.fastq.gz` (the default for NCBI SRA data), the following code would append your sample names directly to a new config file you could use entitled `config/myDataset/myDataset.yaml`:

```
#make new directory, and copy the template into it
mkdir config/myDataset
cp config/config-template.yaml config/myDataset/myDataset.yaml

#Append sample names to template config, assuming a suffix of `1.fastq.gz` :
for file in `ls intermediate/compute-workflow-intermediate/00-fastq | grep 1.fastq.gz | cut -f1 -d_` ; do printf "  $file : $file\n" ; done >> config/myDataset/myDataset.yaml
```

Now, you need to edit your config file (e.g. `config/myDataset/myDataset.yaml`) to include a unique name for your study (which will be appended to output files), the paths to the databases set up above, the path to `uclust`, and the suffixes for your input files (i.e. what should be stripped off to get the sample identifier). After opening the file in your favourite editor, look for and edit the following lines:

```
suffixR1: "_1.fastq.gz" #NCBI format
suffixR2: "_2.fastq.gz" #NCBI format
cutoff: 0.01 #cutoff for classify pipeline, probably can keep the same
phyloFlashDB: "/home/jesse/databases/phyloFlash/138.1/" #location of phyloFlash database, download with phyloFlash's built-in script
bbsplitDBpath: "/home/db/bbsplit-db/" #Download here: https://osf.io/e65rs/
uclustpath: "/home/jesse/MGPrimerEval-tutorial/bin"
study: <your study here>
```

For example, if you used the above commands to install the databases, the paths would be:

```
/home/<your username>/databases/phyloFlash-db/138.1/
/home/<your username>/databases/bbsplit-db/
```

If you don't need or want to test all the primers specified in the config, just comment them out. If you want to add new primers, you have to provide their location on the 4 different SSU rRNA references (see [this repository](https://github.com/jcmcnch/primer-regions.alignments) for an example of how to do so). Or you can just open a github issue and I can add them to the template config. *NB: If you're adding new primers, don't forget to reverse complement the reverse primer sequence.*

## Running the *Compute* workflow

The simplest invocation would be as follows:

```
conda activate snakemake-env
snakemake --cores <# of cores> --use-conda --snakefile Snakefile-compute.smk --configfile config/myDataset/myDataset.yaml
```

To make sure things are set up properly you can:

- Append `-np` to just print out the commands to be run (good way to check if all the input files can be found)
- Append `--conda-create-envs-only` to just install the necessary software using conda

Once the compute pipeline has completed, I recommend you take a look at the output to make sure you're getting good data back. For example, you can scroll through the output on the command line as follows:

```
#concatenate output and pipe to less to look at the data
cat intermediate/compute-workflow-intermediate/09-summary/* | less
```

You should look at the 6th column to check for the number of QC'd sequences recovered per sample/direction/group/primer. It makes sense to start with  **BACT-NON-CYANO** since this category will typically be the most abundant. It's normal to see zeroes or very low numbers for some categories. For example *Archaea* are just not that abundant in many surface ocean samples. If the data is *all zeroes*, it likely means that `uclust` is not set up properly. You can confirm this by looking into the alignment files as follows:

```
#view alignment files to make sure they're not empty:
cat intermediate/compute-workflow-intermediate/05-pyNAST-aligned/* | less
```

If you want, you can just stop here. But if you want to get more information about which taxa are missed by your primers, you can run the *Classify* workflow (next section).

## Running the *Classify* workflow:

The classify workflow summarizes all the data from the compute workflow. For example, if you had 500 samples it would give you information on the overall patterns summed across all 500 samples.

I recommend first running a small portion of the workflow to generate a summary plot from the *Compute* results:

```
#Activate snakemake and run only the summary plot step
conda activate snakemake-env
#The until flag only runs the pipeline up to a specific step
snakemake --cores <# of cores> --use-conda --snakefile Snakefile-classify.smk --configfile config/tutorial/config.yaml --until plot_compute_results 
```

This will produce boxplots summarizing the percent of environmental sequences perfectly matching individual primers in your environment at 0, 1, and 2-mismatch thresholds. If you have a small amount of input data, or your environment is without one of the groups, it's normal for some of the plots to be blank.

Once you verify (at least some of) the plots look good, you can run the rest of the pipeline as follows:

```
#Activate snakemake and run full pipeline
conda activate snakemake-env
snakemake --cores <# of cores> --use-conda --snakefile Snakefile-classify.smk --configfile config/tutorial/config.yaml
```

## Example data (e.g. if you just want to test/verify the functionality of the pipeline on your system)

The following instructions 

**These steps were tested on a remote server running Ubuntu 16.04 in January 2021.**

I'm going to assume you're familiar with [basic bash command line syntax](https://astrobiomike.github.io/unix/unix-intro), have github installed, and you're using something like `screen` or `tmux` to keep a persistent session alive. I'll also assume you've followed the [snakemake install instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) with the small difference that my conda environment for snakemake is `snakemake-env` not `snakemake`.

### Getting setup: Downloading databases

Our pipeline depends on open-source software, and leans heavily in particular on the wonderful [phyloFlash package](https://github.com/HRGV/phyloFlash/blob/master/README.md) which uses a curated version of the SILVA database to retrieve SSU rRNA from meta'omics data (both 16S and 18S). Before running the tutorial, you will need to download their database. Choose a sensible location, install phyloFlash using conda (the following worked for me: `conda config --set channel_priority true && conda create -n pf sortmerna=2.1b phyloflash && conda config --set channel_priority false`), activate the environment (`conda activate pf`), run their download script in the sensible location (e.g. `phyloFlash_makedb.pl --remote`; this step will take up to a few hours as phyloFlash runs quality-control steps on the database but this only needs to be run once), and then make a note of where the database is located so you can add it into your configuration file later.

There is also a splitting step for which you will need to download the sequence files found [here](https://osf.io/e65rs/). To transform these sequence files into splitting databases for the pipeline, you will need `bbtools` installed ([install with conda here](https://anaconda.org/bioconda/bbmap)). Use the commands in the shell script found in the repository. You don't need to have a `bbmap-env` conda environment - as long as your shell can access `bbtools` when you run the commands you should be fine. You will know the databases have been constructed correctly if the folder now contains three subfolders, which are the databases `bbsplit` will use to partition the SSU rRNA data into the 4 different organismal categories (*Archaea*, *Bacteria*, *Cyanobacteria* + Plastid 16S, *Eukarya*). Once this is finished, make a note of the location so you can put that information in your configuration file.

### Running the workflow

I usually make a new folder for each new dataset I'm analyzing to keep things organized. For the purposes of this tutorial, let's download the repo into a folder called `MGPrimerEval-tutorial` as follows:

`git clone https://github.com/jcmcnch/MGPrimerEval.git MGPrimerEval-tutorial`

1. Now, `cd` into the folder. If you have data to play with, add it into the folder `intermediate/compute-workflow-intermediate/00-fastq/`. If not, clone [this extremely helpful repo](https://github.com/wwood/ena-fast-download) into the tutorial folder (you should now see a `ena-fast-download` subfolder), install [ascp](https://download.asperasoft.com/download/docs/ascp/3.5.2/html/index.html) (which is proprietary software from IBM that downloads things *very fast*) and use the script `./tutorial/download-BGT.sh` to download a small subset of the [BioGEOTRACES metagenomes](https://www.nature.com/articles/sdata2018176). The script will also put the files into the right spot. *Keep in mind, this is still a fair bit of data! If possible, do it on a work server, not your home network unless you have unlimited bandwidth.*
2. The next step is to set up your configuration file for snakemake to read so it knows what samples to analyze. If you're using the BioGEOTRACES data, there is already a configuration file in `config/tutorial/config.yaml`. If not, you can use the template at `config/config-template.yaml`. *Make sure to modify it to point to database locations and suit your needs*. For the tutorial, you only need to provide the locations of the datbases (i.e. edit the lines starting with `phyloFlashDB` and `bbsplitDBpath`). If you are using your own data and the template, you need to add a study name and put your sample IDs at the end, as well as providing the appropriate suffixes for your raw data (i.e. what should be stripped off your raw fastq files to get the sample IDs). You can optionally add/remove primers as you see fit (to remove existing primers, just comment them out from the config). If you want to add new primers, you have to provide their location on the 4 different SSU rRNA references (see [this repository](https://github.com/jcmcnch/primer-regions.alignments) for an example of how to do so).
3. Before you run snakemake, you'll have to activate the conda environment i.e. `conda activate snakemake-env`.
4. Now, run the compute step (upon which the other two modules depend) as follows if you downloaded the tutorial data: `snakemake --cores <# of cores> --use-conda --snakefile Snakefile-compute.smk --configfile config/tutorial/config.yaml`. Snakemake will now automagically install all software dependencies and should begin cranking out data. If you want to test whether everything is working before initiating the run, you can append `-np` to the above command. This will do a dry run and print out some useful output about the steps snakemake plans to execute.

Known issues:

* If you are running into issues with DAG generation (read [snakemake documentation](https://snakemake.readthedocs.io/en/stable/) if you want to know what a DAG is) taking a long time, especially if you have a *lot* of samples, you might need to subset your workflow. You can do this manually (some examples of how to do so are found in the `runscripts` folder), or use a [new batch mode](https://snakemake.readthedocs.io/en/stable/executing/cli.html#dealing-with-very-large-workflows) built into the latest versions of snakemake (not implemented in this workflow).
* bbmap/bbsplit steps sometimes can hang under situations of high RAM use - it will just get stuck at a particular step and not proceed further. To resolve this, just kill (i.e. CTRL-C) and restart your workflow.

A visual demonstration of the compute module:
![Rule graph for initial (Snakefile-compute.smk) steps](https://github.com/jcmcnch/MGPrimerEval/blob/master/images/Snakefile-compute.svg)
