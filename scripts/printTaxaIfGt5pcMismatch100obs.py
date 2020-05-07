#!/usr/bin/env python3

import sys
import csv

#get taxa that have enough sequences to evaluate their mismatches across the dataset
#require at least 50 sequences that miss, and 5% mismatch in order to print
#this will avoid printing sequence-based errors or very rare taxa

for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):
	
	#skip header
	if astrLine[0] != "taxonomic_group":
		
		if float(astrLine[3]) > 50 and float(astrLine[6]) >= 0.05:

			print(astrLine[0])
