#!/usr/bin/env python3.7

from Bio import SeqIO

import argparse

parser = argparse.ArgumentParser(description='This script filters a pyNAST alignment for a region of interest (ROI) specified by the user (i.e. a primer region in E. coli coordinates). It will output a fasta file containing only those records that do not contain gaps or ambiguous bases in a specific region of interest in a multiple sequence alignment. This will only work for pyNAST-like formats, which align to a fixed reference sequence.')

parser.add_argument('--input', help='Your multiple sequence alignment (MSA) in gapped fasta format as produced by PyNAST.')

parser.add_argument('--output', help='An output name for the cleaned MSA with only reads covering the region of interest.')

parser.add_argument('--start', help='Start of the region of interest in E. coli 16S rRNA coordinates')

parser.add_argument('--end', help='End of the region of interest in E. coli 16S rRNA coordinates')

parser.add_argument('--padding',help='Number of padding bases required on edges of ROI before and after the region you specify so that the region is not right at the edge of the read. Default: 5 bases.',default='5')

args = parser.parse_args()

hashCleanSeqs = {}

iStart = int(args.start) - int(args.padding)
iEnd = int(args.end) + int(args.padding)

for record in SeqIO.parse(args.input, "fasta"):

	sequence = str(record.seq).upper()

	if ("-" not in sequence[iStart:iEnd]) and ("N" not in sequence[iStart:iEnd]):
	
		hashCleanSeqs[record.description] = sequence
	
with open(args.output, "w+") as output_file:

	for key, value in hashCleanSeqs.items():

		output_file.write(">" + key + "\n" + value + "\n")
