#!/bin/bash

for item in 04-sorted/*EUK.fa; do 

	filestem=`basename $item .EUK.fa`
	eukCount=`grep -c ">" $item`
	prokCount=`grep -c ">" 04-sorted/$filestem.PROK.fa` 
	totalSeqs=$(python -c "print($eukCount + $prokCount)") 
	eukFrac=`bc <<< "scale=8; $eukCount/$totalSeqs"` 
	printf "$filestem\t$eukFrac\n"

done
