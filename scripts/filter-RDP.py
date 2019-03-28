#!/usr/bin/python

import csv
import sys

for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):

    if float(astrLine[2]) >= 0.5:

        print(astrLine[1])
