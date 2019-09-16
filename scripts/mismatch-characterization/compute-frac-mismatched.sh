#!/bin/bash

targetFileInput=$1 # use snakemake to iterate over target files (which contain taxonomy labels that meet a particular user-defined abundance threshold) generated by previous step, rerunning this script multiple times
#Targets include those taxonomy labels meeting the threshold in both the matching and mismatching sequences

totalFilteredSeqs=`cut -f2 $2` #get of number of total hits at particular mismatch threshold

normalizedHitsCounts=$3

missesCountFile=$4

#print summary header
printf "taxonomic_group\ttaxon_primer_hits\ttaxon_primer_misses\ttaxon_total_observed\ttotal_filtered_seqs\tfraction_total_filtered_sequences_mismatched\tfraction_taxon_mismatched\n"

#Looping through the taxonomic groups defined in the target file
while read line ; do

	#for a particular taxonomic group, get the estimated number of counts found in the original dataset (necessary to normalize since we subsample during taxonomic annotation)
	misses=`awk -v line="$line" '$2 == line {print $1}' $normalizedHitsCounts`

	#Get the number of all mismatched sequences for the same taxonomic group (not normalized; the pipeline takes all the mismatches to maximize the information we can retrieve on the primer mismatches)
	hits=`awk -v line="$line" '$1 == line {print $2}' $missesCountFile`

	#set value to zero if no misses (perfectly matched)
	if [ -z "$misses" ] ; then misses=0 ; fi

	#get total
	total=$(python -c "print($hits + $misses)")

	#if statement to avoid division by zero
	if (( $(echo "$misses > 0" | bc -l) )) ; then

		#Calculate fraction of misses relative to total number of sequences for the taxonomic group in question
		#Essentially trying to answer this question: How quantitatively important is this mismatch? For example, a SAR11 mismatch is going to be very important quantitatively because this is an abundant group.
		missesAsFracOfTotal=`bc <<< "scale=8; $misses/$totalFilteredSeqs"`
	else
		missesAsFracOfTotal=0
	fi

	#if statement to avoid division by zero
	if (( $(echo "$total > 0" | bc -l) )) ; then

		#Calculate fraction of this taxonomic group that is missed by primers
		#What fraction of the sequences from this group are missed?
		fracMismatched=`bc <<< "scale=4; $misses/$total"`
	else
		fracMismatched=0
	fi

	printf "$line\t$hits\t$misses\t$total\t$totalFilteredSeqs\t$missesAsFracOfTotal\t$fracMismatched\n"

done < $targetFileInput
