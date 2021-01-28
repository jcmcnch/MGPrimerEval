#!/usr/bin/env python3.7

from Bio import SeqIO
import csv
import sys

aHitIDs = []

for astrLine in csv.reader(open(snakemake.input[0]), csv.excel_tab):

	aHitIDs.append(astrLine[0].strip())

setHitIDs = set(aHitIDs)

aAllIDs = []

for record in SeqIO.parse(snakemake.input[1], "fasta"):

    aAllIDs.append(record.id.strip())

setAllIDs = set(aAllIDs)

setMissing = setAllIDs - setHitIDs #Those that are missing in the blast results

setPresent = setAllIDs - setMissing

aOutputMissing = []

aOutputPresent = []

for record in SeqIO.parse(snakemake.input[1], "fasta"):

	if record.id in setMissing:

		aOutputMissing.append(record)

for record in SeqIO.parse(snakemake.input[1], "fasta"):

	if record.id in setPresent:

		aOutputPresent.append(record)

SeqIO.write(aOutputMissing, snakemake.output[0], "fasta")

SeqIO.write(aOutputPresent, snakemake.output[1], "fasta")
