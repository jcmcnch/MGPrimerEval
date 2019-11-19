#!/usr/bin/env python

"""
Code taken from user "Jon" on biostars, see: https://www.biostars.org/p/163928/
Also:
https://github.com/nextgenusfs/amptk
"""

import re
import csv
import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='This script checks to see if the primer region plus 5 bases on either side has a minimum qualtiy of 30 or more. If so, it outputs a filtered fastq file.')

parser.add_argument('--info', help='The tab-separated info file from cutadapt. This is where we get the primer region info from.')

parser.add_argument('--fastq', help='Your input fastq file to be filtered by the quality scores of the matching region (i.e. primer).')

args = parser.parse_args()

hashQualCoordinates = {}

#create a dictionary with the relevant quality information from the info file
for astrLine in csv.reader(open(args.info), csv.excel_tab):

	if astrLine[1] != "-1":

		hashQualCoordinates[astrLine[0].strip()] = [astrLine[2], astrLine[3]]

#use the information determined above to slice the fastq file to the primer region
for record in SeqIO.parse(args.fastq, "fastq"):

	if str(record.description) in hashQualCoordinates:

		iStart = int(hashQualCoordinates[record.description][0]) - 5

		iEnd = int(hashQualCoordinates[record.description][1]) + 5

		if (iStart >= 0) and (iEnd <= len(record.seq) - 1): #Account for those that fall outside of padding region

			if min(record.letter_annotations["phred_quality"][iStart:iEnd]) >= 30:

				sys.stdout.write(record.format('fastq'))
