Just want to run the pipeline? [Jump to usage instructions.](https://github.com/jcmcnch/MGPrimerEval#running-the-pipeline-with-your-own-data)

## Table of contents:

1. [Preamble](https://github.com/jcmcnch/MGPrimerEval#1-preamble)
2. [Pipeline Architecture](https://github.com/jcmcnch/MGPrimerEval#2-pipeline-architecture)
3. Detailed Overview:  
  3.1 [Motivating Scientific Questions](https://github.com/jcmcnch/MGPrimerEval#31-motivating-scientific-questions)  
  3.2 [Input Requirements](https://github.com/jcmcnch/MGPrimerEval#32-input-requirements)  
  3.3 [Recommended Operating Systems](https://github.com/jcmcnch/MGPrimerEval#33-recommended-operating-systems)  
  3.4 [Overview of Pipeline Steps](https://github.com/jcmcnch/MGPrimerEval#34-overview-of-pipeline-steps)  
  3.5 [Expected Output Files](https://github.com/jcmcnch/MGPrimerEval#35-expected-output-files)  
  3.6 [Open Source Software and Database Dependencies](https://github.com/jcmcnch/MGPrimerEval#36-open-source-software-and-database-dependencies)
4. Running the pipeline with your own data:  
  4.1 [Setting up and activating a snakemake conda environment](https://github.com/jcmcnch/MGPrimerEval#41-setting-up-and-activating-a-snakemake-conda-environment)  
  4.2 [Cloning the repository and adding raw data](https://github.com/jcmcnch/MGPrimerEval#42-cloning-the-repository-and-adding-raw-data)  
  4.3 [Downloading databases](https://github.com/jcmcnch/MGPrimerEval#43-downloading-databases-for-phyloflash-ssu-rrna-splitting-ssu-rrna-classification-with-vsearch-and-adding-uclust-to-your-path)  
  4.4 [Setting up your configuration file](https://github.com/jcmcnch/MGPrimerEval#44-setting-up-your-configuration-file)
5. [Running the Compute workflow](https://github.com/jcmcnch/MGPrimerEval#5-running-the-compute-workflow)
6. [Running the Classify workflow](https://github.com/jcmcnch/MGPrimerEval#6-running-the-classify-workflow)
7. [Running the Compare workflow](https://github.com/jcmcnch/MGPrimerEval#7-running-the-compare-workflow)
8. [Downloading example data](https://github.com/jcmcnch/MGPrimerEval#8-downloading-example-data-eg-if-you-just-want-to-testverify-the-functionality-of-the-pipeline-on-your-system)
9. [Known Issues](https://github.com/jcmcnch/MGPrimerEval#9-known-issues)

## 1. Preamble
Accurate design of oligonucleotide primers for small subunit ribosomal RNA (SSU rRNA) polymerase chain reaction (PCR) amplicon sequencing (or indeed any PCR-based analysis) determines how quantitiative the resulting data is. So far, primers have been mainly designed based on comprehensive and highly-curated reference databases such as SILVA. This has provided important insights into theoretical primer performance and corrected many flaws. However in past primer evaluations, the reliance on full-length references and giving each sequence equal weight can lead to a distorted perspective on the actual extent of matches and mismatches expected in real samples. Prior approaches did not take into account: 

+ The highly unequal abundances in nature, so for example, a very abundant mismatched taxon has the same weight as a well-matched one occurring at vanishingly small abundance, hence the “% match” to the database could significantly over- or under-estimate the actual “% match” to the field even if all taxa are in the database.
+ Some environments may have abundant taxa poorly represented or unrepresented in these curated reference databases. 

Our new approach as implemented with this pipeline is an attempt to provide an automated, less biased way of evaluating primer performance based on meta-'omics datasets from the environment of interest, avoiding both sources of bias mentioned above.

## 2. Pipeline Architecture

The pipeline is divided into three modules:

* *Compute* workflow: Calculate, from raw (unassembled) metagenome/transcriptome reads, the proportion of reads that perfectly match your primer(s) of interest.
* *Classify* workflow: Using output from the *Compute* workflow, determine which taxa are matched, and which taxa are mismatched. Provide information on primer variants for specific groupings that can allow an investigator to correct biases as appropriate.
* *Compare* workflow: Using output from the *Compute* workflow, compare relative abundances of amplicon sequence variants with those from the *same region* of the 16S/18S molecule. This allows an "apples to apples" comparison between amplicon/metagenomic methods, and can show how well the two methods correspond with one another. 

We have described the results of this analysis already for oceanic ecosystems (see preprint [here](https://www.biorxiv.org/content/10.1101/2020.11.09.375543v1)), but we should note that *our pipeline is agnostic to primer/dataset* and thus could be used on any meta-'omics dataset. We hope the instructions below are enough to get you started. If there are any bugs, or questions, please open a github issue above and we'll do our best to help.

## 3. Detailed Overview

### 3.1 Motivating Scientific Questions

This pipeline is designed to address several related questions at different levels of detail:
+ What fraction of environmental SSU rRNA fragments match oligonucleotide primer sequences\* in a given environment at 0, 1, and 2-mismatch thresholds? *i.e., how well do primers theoretically perform for a given environment / dataset?*
+ What taxa are perfectly matched by my primers? What taxa are not, and therefore likely to be inaccurately quantified using PCR amplicon barcoding methods? *i.e. what are the taxonomic blindspots of a given oligo / primer set?*
+ How can I improve a given oligonucleotide primer to improve its performance on a given dataset/environment? *i.e. can I create an "optimal" primer set for my environment?*
+ How well do results from PCR primers compare with metagenomic reads? *i.e. are the two methods consistent or inconsistent with one another?*

\*Note: the pipeline *does* handle degenerate primers as a query, and returns a perfect match if one of the variants specified in your degenerate sequence matches the target.

### 3.2 Input Requirements:

The only thing you need are raw, *unassembled* paired-end meta'omics data such as metagenomes or metatranscriptomes. They should be compressed in gzip format (suffix = `gz`). Merged read pairs or single-end reads are not currently supported. It is important that data have not been filtered or assembled, since the goal of this pipeline is to recover the underlying environmental patterns with respect to primer matches/mismatches (which assumes that your metagenome/-transcriptome is an accurate representation of the environment in question). 

### 3.3 Recommended operating systems

The recommended operating system to run the pipeline is Linux (tested on Ubuntu 16.04/20.04, and CentOS). It may work on Mac/Windows so long as you can install snakemake and conda, but has not been tested on these systems.

### 3.4 Overview of Pipeline Steps:

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

### 3.5 Expected output files

*NB: By default, the pipeline keeps all intermediate files except for the fastp processed raw reads, but you can change this behaviour by putting `temp()` around any output files you wish to discard. That being said, the processed data files should be considerably smaller than your raw data. I also have a cleanup script in the repository you can use to compress and remove some unnecessary intermediates if you're running out of space (`scripts/compress-cleanup-MGPrimerEval.sh`). This script should only be run after you finish running all modules.*

If you run just the *Compute* workflow, you will get:
- phyloFlash summaries of taxa present in your metagenome/-transcriptome (can be useful to make sure your sample is what you think it is)
- QC'd SSU rRNA fragments
- Sorted SSU rRNA fragments (into categories noted above)
- Aligned SSU rRNA fragments
- Fragments subsetted from alignments to primer regions, and sorted into matching and non-matching at 0, 1, and 2-mismatch thresholds
- Summaries of these data in a tab-separated format suitable for importing into plotting software

If you additionally run the *Classify* workflow using output from above, you will get:
- A graphical summary of the overall results from the *Compute* workflow
- Taxonomic assignments for above matching / mismatching fragments (at the order level)
- Graphical summaries of the taxonomic matches / mismatches
- Tabular summaries across your whole dataset, indicating the proportions of order-level taxonomic groups are matched or mismatched for each primer/group/mismatch threshold (e.g. summarizing mismatches to *Archaea* for the primer 515Y at a 0-mismatch threshold)
- Tabular summaries of the primer variants, both for the groups indicated above *and* for taxa that have abundant mismatches (so you can correct these mismatches or at least know what they are, even if the taxa are rare)

The *Compare* workflow (only if you have paired metagenomes and amplicon sequences you want to intercompare):
- SSU rRNA fragments subsetted to the primer region
- A BLASTn-based comparison between MG SSU rRNA fragments and amplicon sequence variants (using ASVs as a BLAST database and the MG SSU rRNA as query)
- A direct intercomparison between taxonomic groups found in MG SSU rRNA and ASVs *from the same sample*, summarized in graphical and tabular format (includes R^2 values of relative abundances; see manuscript text for more details)

### 3.6 Open-source software and database dependencies

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

## 4. Running the pipeline with your own data

The following are instructions to get the pipeline set up for your own datasets. There are also [instructions below for downloading and processing example data](https://github.com/jcmcnch/MGPrimerEval#8-downloading-example-data-eg-if-you-just-want-to-testverify-the-functionality-of-the-pipeline-on-your-system) if you just want to test the mechanics and make sure it runs on your system. You will still need to follow most of the setup below, with the exception of adding your sample names to the config (a config is already provided with the sample names for these example data).

These instructions assume you have familiarity with [basic bash command line syntax](https://astrobiomike.github.io/unix/unix-intro), have github installed, and you're using something like `screen` or `tmux` to keep a persistent session alive for remote servers. 

### 4.1 Setting up and activating a snakemake conda environment

Assuming you have the [python3 version of miniconda](https://conda.io/en/latest/miniconda.html) installed, install snakemake into its own environment and activate it as follows:

```
#install mamba, which is faster than conda
conda install -c conda-forge mamba

#create an environment named snakemake-env
mamba create -c conda-forge -c bioconda -n snakemake-env snakemake
conda activate snakemake-env
```

([Source for install instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html))

### 4.2 Cloning the repository and adding raw data

Now, clone the repo into a new folder we'll call `myDataset` and enter that folder.

```
git clone https://github.com/jcmcnch/MGPrimerEval.git myDataset
cd myDataset
```

Link your raw data into the input folder (the `ln -s` "softlink" prevents data duplication). For example:

`ln -s /full/path/to/your/data/*gz intermediate/compute-workflow/00-fastq/`

### 4.3 Downloading databases for phyloFlash, SSU rRNA splitting, SSU rRNA classification with VSEARCH, and adding `uclust` to your path

**4.3.1 PhyloFlash Database for Retrieving SSU rRNA fragments**

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

**4.3.2 Database for Splitting SSU rRNA fragments**

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

**4.3.3 Databases for classifying SSU rRNA fragments with VSEARCH**

**NB: If these databases are not set up correctly, the *Classify* workflow will run but will generate empty output files.** 

The following code will download the `udb`-formatted files you can use for the *Classify* workflow (not needed if you only want to run the *Compute* workflow):

```
#make directory, enter it
mkdir -p ~/databases/VSEARCH_db/ ; cd ~/databases/VSEARCH_db

#download database files from OSF
for item in 25a8b znrv8 ; do curl -O -J -L https://osf.io/$item/download ; done

#get full paths for config file
ls $PWD/*udb
```

**4.3.4 Getting the `uclust` executable:**

**Please note that the pipeline will still run if `uclust` is not set up correctly, but will produce empty output files for the alignment step (meaning you won't get any results). So double-check that the path you provide in the config below is accurate.**

The alignment steps in this pipeline currently depend on `pyNAST`, which also depends on `uclust`. However, the `uclust` executable is not available through standard repositories as it is not open-source. You may have access to `uclust` (e.g. from an older install of `qiime`), but you can also just email me at mcnichol at alum dot mit dot edu and I'll send you the binary. I have [been given permission](https://github.com/biocore/pynast/issues/21) to distribute the executable by email by the author.

All you need to do is put the binary in a sensible location, and make a note of the full path to add to the config file (next section). This way, when you run the workflow, this executable will be found. For example:

```
#enter your analysis folder
cd myDataset

#make a directory for binaries, enter directory
mkdir bin ; cd bin

#copy uclust into this directory
cp /location/of/uclust .

#make executable
chmod a+x uclust

#get full path to put into the config file (below)
pwd
```

### 4.4 Setting up your configuration file

The template configuration file comes pre-set with a number of primers that we tested in our study. If you just want to test these primers on your samples, all you have to do is add your sample names at the end in the format `  sample : sample`. I suggest making a new folder and config file for your analysis to keep things organized. If your forward reads end with `_1.fastq.gz` (the default for NCBI SRA data), the following code would append your sample names directly to a new config file you could use entitled `config/myDataset/myDataset.yaml`:

```
#make new directory, and copy the template into it
mkdir config/myDataset
cp config/config-template.yaml config/myDataset/myDataset.yaml

#Append sample names to template config, assuming a suffix of `1.fastq.gz` :
for file in `ls intermediate/compute-workflow/00-fastq | grep 1.fastq.gz | cut -f1 -d_` ; do printf "  $file : $file\n" ; done >> config/myDataset/myDataset.yaml
```

Now, you need to edit your config file (e.g. `config/myDataset/myDataset.yaml`) to include a unique name for your study (which will be appended to output files), the paths to the databases set up above, the path to `uclust`, and the suffixes for your input files (i.e. what should be stripped off to get the sample identifier). After opening the file in your favourite editor, look for and edit the following lines:

```
suffixR1: "_1.fastq.gz" #NCBI format
suffixR2: "_2.fastq.gz" #NCBI format
cutoff: 0.01 #cutoff for classify pipeline, probably can keep the same
phyloFlashDB: "/home/jesse/databases/phyloFlash/138.1/" #location of phyloFlash database, download with phyloFlash's built-in script
bbsplitDBpath: "/home/db/bbsplit-db/" #Download here: https://osf.io/e65rs/
uclustpath: "/home/jesse/MGPrimerEval-tutorial/bin"
VSEARCHudbPath: "/home/db/VSEARCH/silva132_99_sintax.udb"
PhytoRefUdbPath: "/home/db/PhytoRef/PhytoRef_plus_Cyano.udb"
study: <your study here>
```

For example, if you used the above commands to install the databases, the paths would be:

```
/home/<your username>/databases/phyloFlash-db/138.1/
/home/<your username>/databases/bbsplit-db/
/home/<your username>/myDataset/bin/
/home/<your username>/databases/VSEARCH_db/silva132_99_sintax.udb
/home/<your username>/databases/VSEARCH_db/PhytoRef_plus_Cyano.udb
```

If you don't need or want to test all the primers specified in the config, just comment them out. If you want to add new primers, you have to provide their location on the 4 different SSU rRNA references (see [this repository](https://github.com/jcmcnch/primer-regions.alignments) for an example of how to do so). Or you can just open a github issue and I can add them to the template config. *NB: If you're adding new primers, don't forget to reverse complement the reverse primer sequence.*

## 5. Running the *Compute* workflow

The simplest invocation would be as follows:

```
conda activate snakemake-env
snakemake --cores <# of cores> --use-conda --snakefile Snakefile-compute.smk --configfile config/myDataset/myDataset.yaml
```

To make sure things are set up properly you can:

- Append `-np` to just print out the commands to be run (good way to check if all the input files can be found)
- Append `--conda-create-envs-only` to just install the necessary software using conda

Once the compute pipeline has completed, I recommend you take a look at the output to make sure you're getting good data back. For example, you can view the output on the command line as follows:

```
#concatenate output and pipe to less to look at the data
cat output/compute-workflow/09-summary/* | less
```

You should look at the 6th column to check for the number of QC'd sequences recovered per sample/direction/group/primer. It makes sense to start with  *BACT-NON-CYANO* since this category will typically be the most abundant. It's normal to see zeroes or very low numbers for some categories. For example *Archaea* are just not that abundant in many surface ocean samples. If the data is *all zeroes*, it likely means that `uclust` is not set up properly. You can confirm this by looking into the alignment files as follows:

```
#view alignment files to make sure they're not empty:
cat intermediate/compute-workflow/05-pyNAST-aligned/* | less

#you should see something like this on your screen if the alignment was successful (will be empty if the alignment failed):
>SRR5720219.765562.1 NS500496_227_HVGG2BGXX:1:11306:17362:3569 length=150 1..68
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GGATGGGCCTGCGGCGTATCAGGTTGTAGGGGGTGTAATGTACCCTCTAGCCTTCGACGCGTACGGGT------------------------------------------------------------------------------------------------------
```

If you want to get more information about which taxa are missed by your primers, you can run the *Classify* workflow (next section). The first step below also produces some graphical output that may be helpful determining whether it's worth proceeding further (for example, if your primers are already nearly perfect for your environment in question, then you may not care to run the next steps).

## 6. Running the *Classify* workflow:

The classify workflow summarizes all the data from the compute workflow. For example, if you had 500 samples it would give you information on the overall patterns summed across all 500 samples. It is mainly useful for finding out which taxa might be primarily responsible for mismatches, and provides information on how to modify primers to fix these issues.

I recommend first running a small portion of the workflow to generate a summary plot from the *Compute* results:

```
#Activate snakemake
conda activate snakemake-env

#The until flag only runs the pipeline up to a specific step
#In this case, it just concatenates the compute results and plots them up
snakemake --cores <# of cores> --use-conda --snakefile Snakefile-classify.smk --configfile config/tutorial/config.yaml --until plot_compute_results 
```

This will produce 4 boxplots (in your base analysis directory) summarizing the percent of environmental sequences perfectly matching individual primers in your environment at 0, 1, and 2-mismatch thresholds. If you have a small amount of input data, or your environment lacks some of the groups, it's normal for some plots to be blank.

Once you verify (at least some of) the plots look good, you can run the rest of the pipeline without the `--until` flag:

```
#Activate snakemake and run full pipeline
conda activate snakemake-env
snakemake --cores <# of cores> --use-conda --snakefile Snakefile-classify.smk --configfile config/tutorial/config.yaml
```

If you want taxon-specific information on primer variants, you need to run the following script:

```
#run a bash script to generate order-level summaries of primer mismatches
#the bash script runs a python script which uses the environment specified in envs/biopython-env.yaml
#the name of the conda environment activated in the bash script is
#automatically generated by snakemake, and may not be the same on your system
./scripts/make-taxa-mismatch-summaries.sh
```

**Note: this script is only generates output for order-level taxa that have greater than 5% total mismatches and 50 total sequences across your dataset. This will therefore only work for deeply-sequenced samples or if you have many replicates with moderate coverage.**

The *Classify* output will be found in `output/classify-workflow/`. If you want to learn how we tackled the primer re-design, you can find detailed information [here](https://osf.io/gr4nc/). 

Here is a brief summary of the contents of the folders:
- `output/classify-workflow/plots/matchVSmismatch-barplots` is a good place to start. It contains plots of the community composiiton of the mismatched and matched reads at 0, 1, and 2-mismatch thresholds. It will give you a general idea which organisms may be over-represented in the mismatches and thus might merit some corrections. Some may be empty if there is insufficient data. 
- `output/classify-workflow/overall-summaries` contains detailed information you may wish to peruse further. For example:
	+ Files ending in `mismatch.summary.tsv` are tab-separated spreadsheets that show statistics about mismatches for order-level taxa. This will give more quantitative information complementary to the plots mentioned above.
	+ Files ending in `taxonFracMismatched.0-2mm.tsv` show the fraction of mismatches across the different mismatch thresholds. This would essentially tell you how hard it will be to improve your primer for each taxon. For example, if your taxon is not well covered at 0-mismatches but well-covered at 1-mismactches then it would only require minor modifications.
	+ Files ending in `aln.summary.tsv` show the relative abundances for the *matching* variants at each threshold and their relative abundances. If your goal is to improve your primer, I suggest you start at the 2-mismatch threshold since it will contain things not matched by your primer. Variants found in the 0-mismatch files are already covered perfectly by your primer.
	+ Files ending in `aln.fasta` contain fasta files of the actual variants identified, and could be useful for plotting/visualization.
- `output/classify-workflow/taxa-mismatch-summaries` contains the output from the bash script mentioned above. It could be useful for making modifications to your primers that are taxa-specific. For example, I used this to make sure that mismatches to specific taxa such as the *Ectothiorhodospirales* were corrected in our new primer design.
- `filtered-0-mismatches`, `normalized-summaries`, `summary-mismatch-overlap-primer-pairs`, and `pasted-summaries` were used for generating the summary figure in the paper (showing coverage as function of primer pair, not individual primers). You can ignore them unless you want to reproduce our figure.

## 7. Running the *Compare* workflow

To run this workflow, you need to provide several files (examples that are in the github repository are shown in parentheses; note that the formatting needs to be exactly the same in order for it to run correctly):

1. An ASV table that has been transformed to relative abundances (e.g. `config/compare/200514_ASV_info/DADA2/PROKs/200519_GA03-GP13_all-16S-seqs.with-tax.proportions.tsv`).
2. Sequences for the ASVs (e.g. `config/compare/200514_ASV_info/DADA2/PROKs/200519_GA03-GP13_all-16S-seqs.with-tax.proportions.fasta`)
3. A tab-separated file that indicates which metagenomic files correspond with which ASV sample (e.g. `config/compare/GA03-GP13-sample-SRA.tsv`)
4. A config file with information on the parameters you want to give the pipeline (e.g. `config/compare/config-tutorial.yaml`)

If you've downloaded the example data (next section), the compare workflow can be run with the following script:

`./runscripts/compare/16s-dada2-97pc-tutorial.sh`

The output can be found, you guessed it, in the `output/compare-workflow` folder. Since this comparison has quite a few "knobs" you can turn, I've made it so that the output folder is named to record parameters specified by the user (so you don't forget later). For example, the output from today's test was named:

`output/compare-workflow/2021-03-02_dada2-old-data_blastnPcID-97_minAbund-0.01_vs_MG/`

Within this folder, you'll find two subfolders:

- `07-MG-vs-ASV-plots` contains plots of the metagenomic relative abundances plotted against the ASV relative abundances on both linear and log scales.
- `07-MG-vs-ASV-stats` contains some statistics about the correlation. Interpret these statistics with caution - they are intended to be a data exploration tool, not an authoritative description. For example, they could help identify samples that differ a lot between MG and ASV-based methods, and could point towards ways to optimize methods.

## 8. Downloading example data (e.g. if you just want to test/verify the functionality of the pipeline on your system)

The following instructions explain how to download 10 files from the [BioGEOTRACES metagenomes](https://www.nature.com/articles/sdata2018176) that can be used with configuration files and template scripts to test the functionality of the pipeline on your system.

1. If you haven't done so, clone the repo into a folder named `myDataset` and enter it:

```
git clone https://github.com/jcmcnch/MGPrimerEval.git myDataset
cd myDataset
```
 
2. Now, clone [this extremely helpful repo](https://github.com/wwood/ena-fast-download) into your analysis folder (you should now see a `ena-fast-download` subfolder).

3. Install [ascp](https://download.asperasoft.com/download/docs/ascp/3.5.2/html/index.html) on your system, which is proprietary software that is a dependency of `ena-fast-download`. You may find [these instructions](https://gist.github.com/mfansler/71f09c8b6c9a95ec4e759a8ffc488be3) helpful for installing `ascp`.

4. Run the script `./tutorial/download-BGT.sh` to download the BioGEOTRACES metagenomes with `ena-fast-download`. *Keep in mind, this is still a fair bit of data! If possible, do it on a work server, not your home network unless you have unlimited bandwidth.*

5. The shell script will put the downloaded files in the proper place (i.e. `intermediate/compute-workflow/00-fastq/`). So once you have set up the databases and configuration file (in this case, make sure to edit `config/tutorial/config.yaml`; [see above for instructions](https://github.com/jcmcnch/MGPrimerEval#downloading-databases-for-phyloflash-ssu-rrna-splitting-ssu-rrna-classification-with-vsearch-and-adding-uclust-to-your-path)), all you need to do to run the *Compute* workflow is invoke the `run_tutorial.sh` script found in the base directory. The remaining steps can be run as noted above, just make sure to substitute your configuration file. 

## 9. Known issues:

* If paths are not set up correctly for `uclust` and the `VSEARCH` databases, the pipeline will still run but produce empty output files.
* As mentioned above, some output plots will be empty if there is insufficient data.
* If you are running the pipeline on a cluster that has a job submission system, you may need to `source ~/.bashrc` before submitting your job to get `conda` to be recognized.
* If you are running into issues with DAG generation taking a long time (read [snakemake documentation](https://snakemake.readthedocs.io/en/stable/) if you want to know what a DAG is), especially if you have a *lot* of samples, you might need to subset your workflow. You can do this manually (some examples of how to do so are found in the `runscripts` folder), or use a [new batch mode](https://snakemake.readthedocs.io/en/stable/executing/cli.html#dealing-with-very-large-workflows) built into the latest versions of snakemake (not implemented in this workflow).
* bbmap/bbsplit steps sometimes can hang under situations of high RAM use - it will just get stuck at a particular step and not proceed further. To resolve this, just kill (i.e. CTRL-C) and restart your workflow.

A visual demonstration of the compute module:
![Rule graph for initial (Snakefile-compute.smk) steps](https://github.com/jcmcnch/MGPrimerEval/blob/master/images/Snakefile-compute.svg)
