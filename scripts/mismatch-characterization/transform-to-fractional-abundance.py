#!/usr/bin/env python3

import csv
import sys

sumTotal = 0

for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):

	sumTotal += int(astrLine[0])

for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):

	fracAbund = float(astrLine[0]) / float(sumTotal)

#	if fracAbund >= float(sys.argv[2]):

	print("\t".join([astrLine[1].strip(), str(fracAbund)]))
