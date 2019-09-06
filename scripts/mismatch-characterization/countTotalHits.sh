#!/bin/bash

for record in `cat primers.txt`; do

	primer=$record

	for line in `cat groups.txt`; do

		group=$line

		totalFilteredSeqs=`cat ../10-checked/$primer/0-mismatch/SRR*$group*filtered.fastq | grep -c "^@"`

		printf "$primer.$group\t$totalFilteredSeqs\n"

	done

done
