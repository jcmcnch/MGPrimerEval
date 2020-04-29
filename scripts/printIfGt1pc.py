#!/usr/bin/env python3

import sys
import csv

#calc total
iTotalMisses=0
for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):
	#skip header
	if astrLine[5] != "fraction_total_filtered_sequences_mismatched":
		iTotalMisses += int(astrLine[2])

for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):
        #skip header
	if astrLine[5] != "fraction_total_filtered_sequences_mismatched":
		if iTotalMisses > 0:
			if float(int(astrLine[2])/iTotalMisses) >= 0.01 and int(astrLine[2]) > 10:
				print('\t'.join([astrLine[0],astrLine[5]]))
