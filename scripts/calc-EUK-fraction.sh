#!/bin/bash

for item in intermediate/compute-workflow-intermediate/04-sorted/*EUK.fastq; do 

	filestem=`basename $item .EUK.fastq`
	eukCount=`grep -c "^@" $item`
	prokCount=`grep -c "^@" intermediate/compute-workflow-intermediate/04-sorted/$filestem.PROK.fastq` 
	totalSeqs=$(python -c "print($eukCount + $prokCount)") 
	eukFrac=`bc <<< "scale=8; $eukCount/$totalSeqs"` 
	printf "$filestem\t$eukFrac\n"

done
