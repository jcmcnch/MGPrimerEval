#!/usr/bin/env python3

import sys
import csv

sumFWDworstCase=0
hashFWD={}
for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):
	hashFWD[astrLine[0]] = astrLine[1]
        #the sum of all the percentages for the taxa passing the threshold
	sumFWDworstCase+=float(astrLine[1])

sumREVworstCase=0
hashREV={}
for astrLine in csv.reader(open(sys.argv[2]), csv.excel_tab):
	hashREV[astrLine[0]] = astrLine[1]
	#the sum of all the percentages for the taxa passing the threshold
	sumREVworstCase+=float(astrLine[1])

aFWDtaxa = []
aREVtaxa = []
setTaxa = set()
for key in hashFWD.keys():
	setTaxa.add(key)
	aFWDtaxa.append(key)
for key in hashREV.keys():
	setTaxa.add(key)
	aREVtaxa.append(key)

for taxon in setTaxa:
	if taxon in aFWDtaxa and taxon in aREVtaxa:
		if hashFWD[taxon] > hashREV[taxon]:
			hashFWD[taxon] = float(hashFWD[taxon]) - float(hashREV[taxon])
		elif hashFWD[taxon] < hashREV[taxon]:
			hashREV[taxon] = float(hashREV[taxon]) - float(hashFWD[taxon])

sumFWDavgCase=0
for value in hashFWD.values():
	sumFWDavgCase += float(value)
	
sumREVavgCase=0
for value in hashREV.values():
        sumREVavgCase += float(value)

if float(sumFWDworstCase + sumREVworstCase) > 0:
	print(float(sumFWDavgCase + sumREVavgCase) / float(sumFWDworstCase + sumREVworstCase) )
else:
	print("NA")
