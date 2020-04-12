#!/bin/bash

mkdir -p filtered

for group in ARCH BACT-CYANO BACT-NON-CYANO EUK ; do
	
	for primer in 515Y 926R ; do

		for study in bioGEOTRACES Guaymas_sediments Malaspina_Acinas_2019 MBARI-bloom Saanich_inlet; do

			for item in tsv-summaries/*$study*$group*$primer*; do

				./02-printIfGt1pc.py $item

			done | sort -k2 -r >> filtered/$study.$group.$primer.gt1pc.gt10obs.txt

		done

	done

done
