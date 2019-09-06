#!/bin/bash

searchString=$1 #a string based off of each unique combo of primer and group

input=$2 #input tsv that has been filtered by abundance

cat $input | cut -f1 | sort | uniq > $3 #generate target file (groups to quantify)

for mismatch in 0-mismatch 1-mismatch 2-mismatch; do

	while read line ; do

		totalFilteredSeqs=`awk -v searchString="$searchString" -F$'\t' 'BEGIN {OFS = FS} $1 == searchString {print $2}' totalFilteredHits.tsv`

		hits=`awk -v line="$line" '$1 == line {print $2}' intermediate/$filestem.0-mismatch.hits.all.order.counts.normalized.tsv`
		misses=`awk -v line="$line" '$2 == line {print $1}' intermediate/$filestem.$mismatch.nohits.all.order.counts.tsv`	
		if [ -z "$misses" ] ; then misses=0 ; fi

		total=$(python -c "print($hits + $misses)")

		if (( $(echo "$misses > 0" | bc -l) )) ; then

                        missesAsFracOfTotal=`bc <<< "scale=8; $misses/$totalFilteredSeqs"`
		else
			missesAsFracOfTotal=0
		fi

	
		if (( $(echo "$total > 0" | bc -l) )) ; then

			fracMismatched=`bc <<< "scale=4; $misses/$total"`
		else
			fracMismatched=0
		fi

		printf "$line\t$hits\t$misses\t$total\t$missesAsFracOfTotal\t$fracMismatched\n" >> output/$filestem.$mismatch.summary.tsv

	done < intermediate/$filestem.targets

done


for item in output/*; do 

	cat $item | sort -h -k2 | sponge $item

done
