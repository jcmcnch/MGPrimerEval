#!/bin/bash

source activate opedia-env #an environment with pandas installed in python3

mkdir mismatches-summary

for item in subsetted-info-files/*; do 

	alignmentname=mismatches-summary/`basename $item .info`.aln.fasta

	summaryname=mismatches-summary/`basename $item .info`.summary.tsv

	../scripts/make-mismatch-alignments-and-summarize.py --info $item --alignmentout $alignmentname --summaryout $summaryname

done
