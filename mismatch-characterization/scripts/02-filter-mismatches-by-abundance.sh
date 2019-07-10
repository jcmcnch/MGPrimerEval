#!/bin/bash

minAbund=0.01

for item in intermediate/*.nohits.all.order.counts.tsv; do 

	outfile=intermediate/`basename $item .tsv`.min$minAbund.tsv

	scripts/filter-by-fractional-abundance.py $item $minAbund > $outfile

done
