#!/bin/bash

mkdir tax-intermediate

for targetfile in intermediate/*targets; do 

	filestem=`basename $targetfile .targets`

	while read line; do

		./scripts/outputMatchingIDs.py --query $line --input intermediate/$filestem.0-mismatch.nohits.all.tax --output tax-intermediate/$filestem.$line.0-mismatch.nohits.all.ids 

	done < $targetfile

done
