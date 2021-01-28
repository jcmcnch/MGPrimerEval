#!/usr/bin/env python3.7

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

parser = argparse.ArgumentParser(description='This script grabs a very small region of interest from a fastq file that corresponds to part of the 16s/18s gene identified using graftM (HMMSEARCH) and aligned to the appropriate SSU reference sequence with pyNAST.')

parser.add_argument('--strand', help='A tab-separated file with fastq identifiers and strand inferred by HMMSEARCH.')

parser.add_argument('--lenROI', help='The length of your region of interest. E.g. the number of bases of the primer.')

parser.add_argument('--padding', help='The number of bases to pad the resulting fastq slice on either end.')

parser.add_argument('--start', help='Start of the region of interest in SSU rRNA coordinates appropriate to the reference')

parser.add_argument('--fasta', help='Your fasta file containing the sequences to be sliced.')

parser.add_argument('--fastq', help='Your fastq file containing the sequences to be sliced.')

args = parser.parse_args()

hashStrand = {}

#create a dictionary with the strand information to later reverse complement fastqs as necessary
for astrLine in csv.reader(open(args.strand), csv.excel_tab):

	hashStrand[astrLine[0].strip()] = astrLine[1].strip()

iROILength = int(args.lenROI)

iPadding = int(args.padding)

iStart = int(args.start) - 1 - iPadding

hashCoordinates = {}

#parse input fasta file to determine the number of actual nucleotides before the start site (i.e. ignoring gaps)
for record in SeqIO.parse(args.fasta, "fasta"):

	strSequenceSlice = str(record.seq).upper()[:iStart]

	iLeadingBases = len(re.findall("[ACTG]", strSequenceSlice))

	hashCoordinates[record.id] = iLeadingBases

#use the information determined above to slice the fastq file to the primer region, taking into account strandedness
for record in SeqIO.parse(args.fastq, "fastq"):

	if record.id in hashCoordinates.keys(): #To account for those culled from the alignment

		iPrimerBeg = hashCoordinates[record.id]

		iPrimerEnd = iPrimerBeg + iROILength + (iPadding * 2)

		if hashStrand[record.id] == "+":

			PrimerRegion = record[iPrimerBeg:iPrimerEnd]

			sys.stdout.write(PrimerRegion.format('fastq'))

		elif hashStrand[record.id] == "-":

			PrimerRegion = record.reverse_complement(id=True, name=True, description=True)[iPrimerBeg:iPrimerEnd]

			sys.stdout.write(PrimerRegion.format('fastq'))
