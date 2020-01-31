#!/usr/bin/env python

from Bio import SeqIO, Seq
from Bio.Alphabet import generic_dna

import argparse

parser = argparse.ArgumentParser(description='This script takes as input a pyNAST alignment and corresponding sequences from a fastq file, checks whether any of the sequences have been reverse complemented during pyNAST alignment, then outputs a new fastq file that has been reverse complemented according to the information from the pyNAST file. This way all the sequences are in the proper orientation for tools such as cutadapt that do not do revcomping when searching for patterns. Script is brittle, breaking if headers in pyNAST alignment do not have 4 elements with the last one being added by pyNAST.')

parser.add_argument('--inpynast', help='Your multiple sequence alignment (MSA) in gapped fasta format as produced by PyNAST.')

parser.add_argument('--infastq', help='The fastq corresponding to the sequences in the pyNAST alignment.')

parser.add_argument('--outfastq', help='An output name for a new fastq that has reads reverse complemented as necessary.')

args = parser.parse_args()

hashSeqs = SeqIO.to_dict(SeqIO.parse(args.inpynast, "fasta"))

open(args.outfastq, "a+") #Makes empty file in case there are no records to parse in the input fastq

for record in SeqIO.parse(args.infastq, "fastq"):

	with open(args.outfastq, "a+") as output_file:

		if record.id in hashSeqs.keys():

			#Parse information added by pyNAST - the element may differ depending on your header formatting
			#Comes in this format "RC:$start..$end", where RC: is not always present and start/end indicate the part of the sequence aligned successfully by pyNAST.
			#RC = reverse complemented
			#CHANGE the number for slicing the array below IF YOUR FORMATTING DIFFERS

			pyNASTnfo = hashSeqs[record.id].description.split()[2] #EBI headers different than NCBI so need to adjust here
			degappedSeq = str(hashSeqs[record.id].seq).upper().replace("-", "").strip()

			if "RC" not in pyNASTnfo: #not reverse complemented

				start = int(pyNASTnfo.split("..")[0]) - 1
				end = int(pyNASTnfo.split("..")[1])

				SeqIO.write(record[start:end], output_file, "fastq")

			elif "RC" in pyNASTnfo: #reverse complemented

				start = int(pyNASTnfo.split(":")[1].split("..")[0]) - 1
				end = int(pyNASTnfo.split(":")[1].split("..")[1])

				SeqIO.write(record.reverse_complement(id=True, name=True, description=True)[start:end], output_file, "fastq")
