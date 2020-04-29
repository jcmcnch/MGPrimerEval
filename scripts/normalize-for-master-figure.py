#!/usr/bin/env python3

#example usage: "./scripts/normalize-for-master-figure.py {input} `cat {params.normFactor}`"
import sys
import csv

#take as input the normalization factor calculated in previous step
#if no normalization factor calculated, set to 1
if sys.argv[2] != "NA":
	normFactor=float(sys.argv[2])
else:
	normFactor = 1

#string name
primerPair=str(sys.argv[3])

#array for printing output
outArray=[]

#open pasted file containing summary data from forward and reverse primer for given primer set
for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):
	#require at least 10 observations for both forward and reverse primer to calculate summary statistics. if less than 10, the script will not print output for that sample.
	if int(astrLine[6]) >= 10 and int(astrLine[15]) >=10:
		#assuming completely additive
		worstCaseMissed = ( 1 - float(astrLine[8]) ) + ( 1 - float(astrLine[17]) )
		#account for overlap, using average dataset normalization factor (group-specific)
		avgCase = 1 - ( worstCaseMissed * normFactor )
		#account for possibility of negative values, set to zero if the sum of the missed fractions > 1
		#e.g. if the forward primer misses 60% and the reverse primer also misses 60%, then it would sum to 1.2
		if avgCase < 0:
			avgCase = 0
		worstCase = 1 - worstCaseMissed
		if worstCase < 0:
			worstCase = 0
		#output handling to stdout
		outArray = astrLine
		outArray.append(primerPair)
		outArray.append(str(avgCase))
		outArray.append(str(worstCase))
		print('\t'.join(outArray))
