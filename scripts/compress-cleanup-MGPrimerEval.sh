#!/bin/bash

#do this one only if you've finished the classify workflow
find compute-workflow-intermediate/08-checked/ -name "*fastq" -or -name "*info" -type f | parallel --max-procs 10 -t gzip 

find compute-workflow-intermediate/02-phyloFlash_sifted/ -name "*fastq" | parallel --max-procs 10 -t gzip
find compute-workflow-intermediate/03-low-complexity-filtered/ -name "*fastq" | parallel --max-procs 10 -t gzip
find compute-workflow-intermediate/04-sorted/ -name "*fastq" | parallel --max-procs 10 -t gzip
find compute-workflow-intermediate/05-pyNAST-aligned/ -name "*fastq" | parallel --max-procs 10 -t gzip
find compute-workflow-intermediate/06-subsetted/ -name "*fasta" | parallel --max-procs 10 -t gzip
find compute-workflow-intermediate/07-subsetted-fastq/ -name "*fastq" -type f | parallel --max-procs 10 -t gzip 

rm compute-workflow-intermediate/05-pyNAST-aligned/*log
rm compute-workflow-intermediate/02-phyloFlash_sifted/phyloFlash-other/*sam
rm compute-workflow-intermediate/02-phyloFlash_sifted/phyloFlash-other/*hitstats
