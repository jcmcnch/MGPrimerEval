#!/bin/bash

for item in intermediate/*.hits.all.tax; do

	filestem=`basename $item .hits.all.tax`

	group=`echo $item | cut -d\. -f2`
	primer=`echo $item | cut -d\. -f3`

	totalSeqs=`find ../07-subsetted/ -type f -name "*${group}*${primer}*" -print0 | xargs -0 cat | grep -c ">"`
	subsampledSeqs=`cat $item | wc -l`
	fracSubsampled=`bc <<< "scale=4; $subsampledSeqs/$totalSeqs"`
	
	#printf "$totalSeqs\t$subsampledSeqs\t$fracSubsampled\n"

	while read line; do

		num=`echo $line | awk '{print $1}'`
		tax=`echo $line | awk '{print $2}'`
		normalizedNum=`bc <<< "scale=4; $num/$fracSubsampled"`
		printf "$tax\t$normalizedNum\n"	

	done < intermediate/$filestem.hits.all.order.counts.tsv > intermediate/$filestem.hits.all.order.counts.normalized.tsv

done
