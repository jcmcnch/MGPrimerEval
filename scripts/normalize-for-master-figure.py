#!/usr/bin/env python3

#"./scripts/normalize-for-master-figure.py {input} `cat {params.normFactor}`"
import sys
import csv

#take as input the normalization factor calculated in other script
if sys.argv[2] != "NA":
	normFactor=float(sys.argv[2])
else:
	normFactor = 1

primerPair=str(sys.argv[3])

outArray=[]
#open pasted file
for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):
	#require at least 10 observations for both forward and reverse primer to consider sample
	if int(astrLine[6]) >= 10 and int(astrLine[15]) >=10:
		#assuming completely additive
		worstCaseMissed = ( 1 - float(astrLine[8]) ) + ( 1 - float(astrLine[17]) )
		#account for overlap, using average dataset normalization factor (group-specific)
		avgCase = 1 - ( worstCaseMissed * normFactor )
		worstCase = 1 - worstCaseMissed
		#output handling
		outArray = astrLine
		outArray.append(primerPair)
		outArray.append(str(avgCase))
		outArray.append(str(worstCase))
		print('\t'.join(outArray))
