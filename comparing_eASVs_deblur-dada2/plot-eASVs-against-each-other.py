#!/usr/bin/env python3

import csv
import sys
import pandas as pd

#########################

def eASV_2_hash(table):
	"""
	Convert eASV table to hash (as output from transform-to-proportions.py with the top line stripped).
	"""
	hashASV, aHeaders = {}, []
	with open(table, newline="") as f:
		reader = csv.reader(f, csv.excel_tab)
		aHeaders = next(reader) #Get headers from first line

	#populate dictionary keys
		for item in aHeaders:
			hashASV[item.strip()] = []

		#populate array with rest
		for row in reader:
			for i in range(0,len(aHeaders)):
				hashASV[aHeaders[i]].append(row[i].strip())

	return hashASV

#########################


#Read in blast results
hashConversionTable = {}

for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):
	try:
		hashConversionTable[astrLine[0]].append(astrLine[1])
	except KeyError:
		hashConversionTable[astrLine[0]] = [astrLine[1]]

#Read in the two OTU/eASV tables
hashDADA2 = eASV_2_hash(sys.argv[2])
hashDeblur = eASV_2_hash(sys.argv[3])

aSamples = list(hashDADA2.keys())[2:]
aSamples.remove("S0273") #zero values in DADA2
aSamples.remove("S0362") #zero values in DADA2


#Loop through conversion table, creating new output
hashComparison = {}

for sample in aSamples:
	hashComparison = {}
	for key, array in hashConversionTable.items():
		floatFracDeblur = floatFracDADA2 = 0
		iLocation = hashDeblur["#OTU ID"].index(key) #Get location in list
		floatFracDeblur = float(hashDeblur[sample][iLocation])
		for item in array:
			iLocation = hashDADA2["#OTU ID"].index(item)
			floatFracDADA2 = float(hashDADA2[sample][iLocation]) + floatFracDADA2
		hashComparison[key] = [floatFracDADA2, floatFracDeblur]

	df=pd.DataFrame.from_dict(hashComparison, orient="index")
	outname = "results/correlations/" + sample + ".tsv"
	df.to_csv(outname, sep="\t", header=False)
