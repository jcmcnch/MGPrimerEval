#!/usr/bin/env python3

import csv
import argparse

parser = argparse.ArgumentParser(description='This script subsets an info file from cutadapt given input of a file with sequence identifiers.')

parser.add_argument('--inputids', help='A file with input ids on each line.')

parser.add_argument('--inputinfo', help='An info file that comes from cutadapt that is to be subsetted.')

parser.add_argument('--output', help='A subsetted info file containing only the IDs specified in the --inputids argument.')

args = parser.parse_args()

aInputIDs = []

outhandle = open(args.output, "w+")

for strLine in open(args.inputids):

	aInputIDs.append(strLine.strip())

for astrLine in csv.reader(open(args.inputinfo), csv.excel_tab):

	if astrLine[0] in aInputIDs:

		outhandle.write('\t'.join(astrLine[0:]) + '\n')
