#!/bin/bash

mkdir -p summaries

for group in ARCH BACT-CYANO BACT-NON-CYANO EUK ; do

	for study in bioGEOTRACES Guaymas_sediments Malaspina_Acinas_2019 MBARI-bloom Saanich_inlet; do

		./03-primer-pair-subtract-overlap.py filtered/$study.$group.515Y.gt1pc.gt10obs.txt filtered/$study.$group.926R.gt1pc.gt10obs.txt > summaries/$study.$group.avgCase.txt

	done

done

for item in summaries/*; do

	contents=`cat $item`
	name=`basename $item | cut -f1,2 -d\.`
	printf "$name\t$contents\n"

done > avgCase-relative-to-worst-case.tsv
