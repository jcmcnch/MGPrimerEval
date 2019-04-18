#!/usr/bin/env python3

import csv
import pandas as pd
#from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='This script parses information from a cutadapt info file and produces alignments and a summary of the different sequence variants. Requires biopython and pandas in your environment (can use opedia-env on kraken.usc.edu).')

parser.add_argument('--info', help='The tab-separated info file from cutadapt.')

parser.add_argument('--alignmentout', help='A MSA-formatted fasta for just the ROI (region of interest).')

parser.add_argument('--summaryout', help='A TSV summary of the variants and their relative abundances.')

args = parser.parse_args()

hashSeqs = {}

#create a dictionary with the relevant headers and ROI sequences
for astrLine in csv.reader(open(args.info), csv.excel_tab):

	if astrLine[1] != "-1": #Ignore false positives

		hashSeqs[astrLine[0].strip()] = astrLine[5].strip()

aSeqs = []

for strSeq in hashSeqs.values(): #fill array with sequences

	aSeqs.append(strSeq)

hashSummary = {}

for seq in set(aSeqs): #count sequences

	hashSummary[seq] = str( round( (aSeqs.count(seq) / len(aSeqs)), 3 ) )

summaryDF = pd.DataFrame.from_dict(hashSummary, orient='index', columns=['relative abundance'])
summaryDF = summaryDF.sort_values("relative abundance", axis=0, ascending=False)
summaryDF.to_csv(args.summaryout, encoding='utf-8', sep="\t", index_label="variant")
