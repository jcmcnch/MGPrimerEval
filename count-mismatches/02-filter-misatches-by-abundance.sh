#!/bin/bash

minAbund=0.01

for item in *.nohits.all.order.counts.tsv; do 

	outfile=`basename $item .tsv`.min$minAbund.tsv

	filter-by-fractional-abundance.py $item $minAbund > $outfile

done
