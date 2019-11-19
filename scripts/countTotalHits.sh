#!/bin/bash

for primer in 926R 806RB 515Y V4F V4R V4RB 341F 785R 27F 1389F 1510R ; do

	for group in ARCH BACT-CYANO BACT-NON-CYANO EUK ; do

		totalFilteredSeqs=`cat 10-checked/$primer/0-mismatch/SRR*$group*filtered.fastq | grep -c "^@"`

		printf "$primer\t$group\t$totalFilteredSeqs\n"

	done

done
