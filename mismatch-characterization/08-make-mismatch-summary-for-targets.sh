#!/bin/bash

mkdir tax-intermediate

for targetfile in intermediate-test/*targets; do 

	filestem=`basename $targetfile .targets`

	while read line; do

		cat intermediate-test/$filestem.0-mismatch.nohits.all.tax | sed -re 's/\([0-9]{1}\.[0-9]{2}\)//g' | grep "$line" | cut -f1 > tax-intermediate/$filestem.$line.0-mismatch.nohits.all.ids 

	done < $targetfile

done
