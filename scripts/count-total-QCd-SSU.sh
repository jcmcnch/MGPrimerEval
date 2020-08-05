#!/bin/bash
#merge fwd and reverse
for item in `ls compute-workflow-intermediate/03-low-complexity-filtered/*gz | cut -f1 -d. | sort | uniq` ; do 

	#get name
	ID=`basename $item`
	#count total number of lines
	lines=`zcat $item* | wc -l`
	#since it's a fastq, divide by 4
	total=`expr $lines / 4`
	#print out as tsv
	printf "$ID\t$total\n"

done
