#!/bin/bash

while read line ; do

	./ena-fast-download/ena-fast-download.py $line --ascp-args '-m 1m -l 100m'
        mv *gz compute-workflow-intermediate/00-fastq/

done < tutorial/accessions/GA03-stn 
