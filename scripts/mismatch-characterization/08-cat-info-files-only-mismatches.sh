#!/bin/bash

mismatches="6-mismatch" #This is the arbitrary limit for mismatches I've kept. Anything above this is considered to be a false positive

mkdir catted-mismatch-info-files

for item in `cat primers.txt`; do

	primer=$item

	for item in `cat groups.txt`; do

		group=$item

		cat ../10-checked/$primer/$mismatches/SRR*${group}*.info | awk -F$'\t' 'BEGIN {OFS = FS} $2 >= '1' {print}' > catted-mismatch-info-files/$primer.$group.$mismatches.mismatches.info #output info file only for mismatches ($2 >= '1')

	done

done

