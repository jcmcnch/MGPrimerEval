#!/bin/bash

mkdir output

for item in intermediate/*.0-mismatch.nohits.all.order.counts.min0.01.tsv; do

	filestem=`basename $item .0-mismatch.nohits.all.order.counts.min0.01.tsv`

        cat intermediate/$filestem*0-mismatch*.min0.01.tsv | cut -f1 | sort | uniq > intermediate/$filestem.targets

	for mismatch in 0-mismatch 1-mismatch 2-mismatch; do

		while read line ; do

			hits=`awk -v line="$line" '$1 == line {print $2}' intermediate/$filestem.0-mismatch.hits.all.order.counts.normalized.tsv`
			#echo $hits
			misses=`awk -v line="$line" '$2 == line {print $1}' intermediate/$filestem.$mismatch.nohits.all.order.counts.tsv`	
			if [ -z "$misses" ] ; then misses=0 ; fi
 
			#echo "$misses"
			#if `awk -v line="$line" '$2 == line' intermediate/$filestem.$mismatch.nohits.all.order.counts.tsv`; then
			#	misses=`awk -v line="$line" '$1 == line {print $2}' intermediate/$filestem.$mismatch.nohits.all.order.counts.tsv`
			#else
			#	misses=0
			#fi
			
			total=$(python -c "print($hits + $misses)")
			#echo "$hits"
			#echo "$misses"
			#echo "$total"
			if (( $(echo "$total > 0" | bc -l) )) ; then
				fracMismatched=`bc <<< "scale=4; $misses/$total"`
			else
				fracMismatched=0
			fi

			printf "$line\t$hits\t$misses\t$total\t$fracMismatched\n" >> output/$filestem.$mismatch.summary.tsv

		done < intermediate/$filestem.targets

	done

done

#for item in output/*; do 

#	cat $item | sort -h -k2 | sponge $item

#done
