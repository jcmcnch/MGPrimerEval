#!/bin/bash

#this won't work unless you have the following repository:
#https://github.com/wwood/ena-fast-download
#AND ascp:
#https://download.asperasoft.com/download/docs/ascp/3.5.2/html/index.html

while read line ; do

	./ena-fast-download/ena-fast-download.py $line --ascp-args '-m 1m -l 100m'
        mv *gz intermediate/compute-workflow-intermediate/00-fastq/

done < tutorial/accessions/GA03-KN204-stn20.sra.ids 
