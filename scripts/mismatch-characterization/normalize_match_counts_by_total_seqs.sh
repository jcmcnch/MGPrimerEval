#!/bin/bash

#For matched sequences, we take a subsample prior to classification to reduce computational time
#But this means that in order to estimate the abundance of a particular taxonomic group in the overall sequences, we need to normalize these numbers
#This is done by finding what fraction of original, QC'd sequences were subsampled and then normalizing the final numbers by this fraction

#Input comes from "Snakemake-classify.smk"
inputPath=$1
groupPrimerPattern=$2
subsampledMatchesExpression=$3
countTable=$4

#use grep to calculate the total number of ALL matched sequences identified, in order that we can normalize the subsampled sequence numbers
#Use filtered fastq sequences to get only those passing QC
#Expression counts all sequence records corresponding to group and primer that hit the primer perfectly (i.e. merges all samples)
totalSeqs=`find ${inputPath} -type f -name "*${groupPrimerPattern}.hit.filtered.fastq" -print0 | xargs -0 cat | grep -c "^@"` #use -print0 option for find and -0 option for xargs to handle whitespace if present

#Count total number of subsampled seqs post classification and deconcatenation - this is per group/primer, for all samples
subsampledSeqs=`find classify-workflow-intermediate/02-matches-subsampled/ -type f -name "*${subsampledMatchesExpression}" -print0 | xargs -0 cat | grep -c "^@"`

#Calculate fraction subsampled; loop to handle corner case of 100% sampling
if (( $(echo "$subsampledSeqs > 0" | bc -l) )) ; then

	fracSubsampled=`bc <<< "scale=4; $subsampledSeqs/$totalSeqs"`
else
	fracSubsampled=1 #not technically true but prevents a divide by zero error below, which occurs because there are some empty files generated from the cutadapt step that get propagated through the pipeline
fi

#Now, iterate through the lines of the count file, and normalize numbers
while read line ; do

	num=`echo $line | awk '{print $1}'` #First column is the unnormalized number

	tax=`echo $line | awk '{print $2}'` #Second column is the taxonomic label

	normalizedNum=`bc <<< "scale=4; $num/$fracSubsampled"` #Normalize counts

	printf "$tax\t$normalizedNum\n" #Print the taxonomy and the normalized number, which represents the number of sequences we would expect in the whole dataset (for this group/primer). This allows for subsequent calculations showing the quantitative importance of mismatches

done < ${countTable} #Take input from stdin
