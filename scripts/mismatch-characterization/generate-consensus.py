#!/usr/bin/env python

#Simple script for determining a consensus primer sequence MGPrimerEval tsv output
#takes three arguments after script name:
#1. Input file
#2. Original degenerate primer sequence
#3. Desired abundance cutoff (to exclude rare variants that may be errors or at least are not quantitatively important)

import csv
import sys

#Conversion table for ambiguous nucleotides, degenerate_consensus from bio.motifs doesn't do the trick
hashIUPACambiguousCodes = {
"M": ["A", "C"],
"R": ["A", "G"],
"W": ["A", "T"],
"S": ["C", "G"],
"Y": ["C", "T"],
"K": ["G", "T"],
"V": ["A", "C", "G"],
"H": ["A", "C", "T"],
"D": ["A", "G", "T"],
"B": ["C", "G", "T"],
"N": ["A", "C", "G", "T"]
}

primerSeq = str(sys.argv[2])
cutoff = float(sys.argv[3])

aSeqs = []

for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab): #iterate through table

    if astrLine[1] != "relative abundance" and float(astrLine[1].strip()) >= cutoff: #abundance cutoff
        
        aSeqs.append(astrLine[0].strip())

aSetSeqs = []

for i in range(0, len(aSeqs[0])): #Primer positions

	aPos = []

	for item in aSeqs:

		aPos.append(item[i])

	aSetPos = set(aPos)

	aSetSeqs.append(list(aSetPos))

strSeq = "" #The degenerate sequence to be added to the existing primer

for array in aSetSeqs:

	if len(array) > 1: #Check if more than one base present

		aQuery4Dict = sorted(array) #Sort to compare array against lookup table

		#Find the key (ambiguous nucleotide) matching the sorted array
		strSeq = strSeq + hashIUPACambiguousCodes.keys()[hashIUPACambiguousCodes.values().index(aQuery4Dict)]

	else:
		#Just store non-ambiguous nucleotide
		strSeq = strSeq + array[0]

aSeqConsensus = []
strConsensus = "" #A consensus between the primer and the variants identified from the MG

#Iterate across primer seq
for i in range(0, len(primerSeq)):

	if primerSeq[i] in hashIUPACambiguousCodes.keys():

		aConsensus = hashIUPACambiguousCodes[primerSeq[i]] + aSetSeqs[i]
		aSeqConsensus.append(sorted(list(set(aConsensus))))

	else:

		aConsensus = list(primerSeq[i]) + aSetSeqs[i]
		aSeqConsensus.append(sorted(list(set(aConsensus))))

for i in range(0, len(aSeqConsensus)):

	aQuery4Dict = aSeqConsensus[i]

	if len(aSeqConsensus[i]) > 1:

		strConsensus = strConsensus + hashIUPACambiguousCodes.keys()[hashIUPACambiguousCodes.values().index(aQuery4Dict)]

	else:

		strConsensus = strConsensus + aSeqConsensus[i][0]

print(strConsensus)
