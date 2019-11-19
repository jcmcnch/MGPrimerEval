#!/usr/bin/env python

from Bio import SeqIO

import argparse

parser = argparse.ArgumentParser(description='Given a sliced pyNAST-like alignment and a corresponding unsliced fastq file, this script outputs a sliced fastq file.')

parser.add_argument('--inputfasta', help='Your multiple sequence alignment (MSA) in gapped fasta format as produced by PyNAST.')

parser.add_argument('--inputfastq', help='The corresponding fastq file.')

parser.add_argument('--output', help='An output name for the cleaned MSA with only reads covering the region of interest.')

args = parser.parse_args()

hashFastq = SeqIO.to_dict(SeqIO.parse(args.inputfastq, "fastq"))

hashFasta = SeqIO.to_dict(SeqIO.parse(args.inputfasta, "fasta"))

#Get coordinates for same region in fastq by string.find()

hashFQcoordinates = {}

for key in hashFasta:

	if key in hashFastq:

		slicedSeq = str(hashFasta[key].seq).upper()
		sequence = str(hashFastq[key].seq).upper()
		start = sequence.find(slicedSeq)
		end = start + len(slicedSeq)

		hashFQcoordinates[key] = [start, end]

with open(args.output, "w+") as output_file:

	for key in hashFQcoordinates.keys():

		iStart = hashFQcoordinates[key][0]
		iEnd = hashFQcoordinates[key][1]

		SeqIO.write(hashFastq[key][iStart:iEnd], output_file, "fastq")
