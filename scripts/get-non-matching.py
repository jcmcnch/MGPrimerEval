#!/usr/bin/env python3

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

aOutput = []

for record in SeqIO.parse(snakemake.input[1], "fasta"):

	if record.id in setMissing:

		aOutput.append(record)

SeqIO.write(aOutput, snakemake.output[0], "fasta")

