#!/usr/bin/env python3

import sys
import csv

headerCheck=True

print('\t'.join(["Taxon","Frac mismatched at 0-mismatches","Frac mismatched at 1-mismatches","Frac mismatched at 2-mismatches"]))

for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):
    if headerCheck:
        headerCheck=False
        continue
    if not headerCheck:
        if float(astrLine[3]) > 10 and float(astrLine[10]) > 10 and float(astrLine[17]) > 10 and float(astrLine[6]) > 0.001: #if at least 10 detected sequences found across all mismatches AND more than 0.1% o f taxon is mismatched @ 0-mm
            print('\t'.join([astrLine[0],astrLine[6],astrLine[13],astrLine[20]]))
