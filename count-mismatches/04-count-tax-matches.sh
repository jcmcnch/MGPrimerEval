#!/bin/bash

for item in *.hits.all.tax; do 

	fileout=`basename $item .tax`.order.counts.tsv

	sed -re 's/\([0-9]{1}\.[0-9]{2}\)//g' $item | cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | tail -f -n +2 | awk '{print $1,"\t",$2}' > $fileout

done
