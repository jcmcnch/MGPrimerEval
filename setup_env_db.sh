#!/bin/bash -i

#if you have just installed a fresh system curl may not be installed
sudo apt-get install curl

currdir=$PWD

#Install miniconda3, saying "yes" to all recommended options
#Comment out if you've already installed on your system
#Note this is for Linux, other system information found here:
#https://conda.io/en/latest/miniconda.html
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh
chmod a+x Miniconda3-py39_4.9.2-Linux-x86_64.sh
./Miniconda3-py39_4.9.2-Linux-x86_64.sh
#conda init bash
rm -f Miniconda3-py39_4.9.2-Linux-x86_64.sh
source ~/.bashrc
conda config --set channel_priority strict

#install mamba, which is faster than conda, -y flag says yes to everything
conda install -c conda-forge mamba -y

#create an environment named snakemake-env
mamba create -c conda-forge -c bioconda -n snakemake-env snakemake -y

#Use mamba to create phyloFlash environment
mamba create -c conda-forge -c bioconda --name pf sortmerna=2.1b phyloflash -y

#If you're getting errors, you may need to run `conda update conda` or do a fresh install of miniconda if updating is not easy (sometimes you get all sorts of incompatibilities which can just be solved by a fresh install)
conda activate pf

#change directory to suit your needs
mkdir -p ~/databases/phyloFlash-db/ ; cd ~/databases/phyloFlash-db/

#run database download/construction script
#will take a few hours while the database is downloaded and QC'd
phyloFlash_makedb.pl --remote

#create bbmap-env
mamba create -c bioconda --name bbmap-env bbmap -y

#activate environment
conda activate bbmap-env

#make directory, enter it
mkdir -p ~/databases/bbsplit-db/ ; cd ~/databases/bbsplit-db

#download database files from OSF
for item in kv3xp eux4r npb2k 4qtev s5j6q 5jmkv eahds ; do curl -O -J -L https://osf.io/$item/download ; done

#make databases
chmod u+x make-dbs-bbsplit.sh ; ./make-dbs-bbsplit.sh

#Setup classification databases
#make directory, enter it
mkdir -p ~/databases/VSEARCH_db/ ; cd ~/databases/VSEARCH_db

#download database files from OSF
for item in 25a8b znrv8 ; do curl -O -J -L https://osf.io/$item/download ; done

cd $currdir
#output information to file database.paths to add to your config file
bbsplitpath=`echo ~/databases/bbsplit-db/` ; printf "bbsplitDBpath: \"${bbsplitpath}\"\n" > database.paths
phyloFlashdb=`echo ~/databases/phyloFlash-db/138.1/` ; printf "phyloFlashDB: \"${phyloFlashdb}\"\n" >> database.paths
dbdir=`echo ~/databases/` ; printf "uclustpath: \"${dbdir}\"\n" >> database.paths
silvaudb=`echo ~/databases/VSEARCH_db/silva132_99_sintax.udb` ; printf "VSEARCHudbPath: \"${silvaudb}\"\n" >> database.paths
phytoudb=`echo ~/databases/VSEARCH_db/PhytoRef_plus_Cyano.udb` ; printf "PhytoRefUdbPath: \"${phytoudb}\"\n" >> database.paths

echo "To run the tutorial, you now need to copy the lines in the file \"database.paths\" into your config found in config/tutorial/config.yaml"

ln -s $PWD/test-input/*fastq.gz $PWD/intermediate/compute-workflow/00-fastq/
