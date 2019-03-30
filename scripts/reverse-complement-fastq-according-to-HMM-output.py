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

parser = argparse.ArgumentParser(description='This script reverse complements individual records from a fastq file based on output from graftM/HMMSEARCH and produces a file that can be passed directly to cutadapt for adapter checking.')

parser.add_argument('--strand', help='A tab-separated file with fastq identifiers and strand inferred by HMMSEARCH.')

parser.add_argument('--fastqinput', help='Your input fastq file to be reverse complemented as appropriate.')

args = parser.parse_args()

hashStrand = {}


#create a dictionary with the strand information to later reverse complement fastqs as necessary
for astrLine in csv.reader(open(args.strand), csv.excel_tab):

	hashStrand[astrLine[0].strip()] = astrLine[1].strip()


#use the information determined above to slice the fastq file to the primer region, taking into account strandedness
for record in SeqIO.parse(args.fastqinput, "fastq"):

		if hashStrand[record.id] == "+":

			sys.stdout.write(record.format('fastq'))

		elif hashStrand[record.id] == "-":

			revcomp = record.reverse_complement(id=True, name=True, description=True)

			sys.stdout.write(revcomp.format('fastq'))
