#!/usr/bin/env python3

import csv
import argparse
import re

parser = argparse.ArgumentParser(description='This script checks for matches of a given taxonomy in a tab-separated file, printing if the taxonomy matches. It uses regular expressions to remove confidence values from the taxonomy string.')

parser.add_argument('--query', help='A query string to search for.')

parser.add_argument('--input', help='The input file containing sequence identifiers in the first tab-separated column and then the taxonomy string in the second.')

parser.add_argument('--output', help='Your output file with the sequence identifiers that match your taxonomy of interest.')

args = parser.parse_args()

outhandle = open(args.output, "w+")

for astrLine in csv.reader(open(args.input), csv.excel_tab):

	strTax = re.sub(r'\(.*?\)', '', astrLine[1]).strip()

	strTaxSubset = ','.join(strTax.split(',')[0:4]) #stripping off (potentially) uninformative taxonomic information

	strQuery = str(args.query).strip()

	if strTaxSubset == strQuery :
		
		outhandle.write(astrLine[0] + "\n")
