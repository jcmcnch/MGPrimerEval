#!/bin/bash

for item in *.hits.all.order.counts.min0.01.tsv; do

	filestem=`basename $item .hits.all.order.counts.min0.01.tsv`

	cat $filestem*.min0.01.tsv | cut -f1 | sort | uniq > $filestem.targets

	while read line ; do

		hits=`grep "$line" $filestem.hits.all.order.counts.normalized.tsv | cut -f2`
		misses=`grep "$line" $filestem.nohits.all.order.counts.tsv | cut -f1`
		total=$(python -c "print $hits + $misses")
		fracMismatched=`bc <<< "scale=4; $misses/$total" || fracMismatched=0`
		printf "$line\t$fracMismatched\n" >> $filestem.summary.tsv

	done < $filestem.targets

done
