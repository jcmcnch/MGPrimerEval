#!/usr/bin/env python3

import csv, sys

hashHits = {}
hashMisses = {}
hitTotal = 0
missTotal = 0
taxaSet = set()

#Conversion table between SRAid and sample ID
header=True
for astrLine in csv.reader(open(sys.argv[1]), csv.excel_tab):
	if header:
		header=False
		continue
	if not header:
		hashHits[astrLine[0]] = astrLine[1]
		hitTotal += float(astrLine[1])

for key in hashHits.keys():
	hashHits[key] = float(hashHits[key]) / hitTotal
	taxaSet.add(key)

header=True
for astrLine in csv.reader(open(sys.argv[2]), csv.excel_tab):
	if header:
		header=False
		continue
	if not header:
	        hashMisses[astrLine[0]] = astrLine[1]
        	missTotal += float(astrLine[1])

for key in hashMisses.keys():
        hashMisses[key] = float(hashMisses[key]) / missTotal
        taxaSet.add(key)

hashOutputTable = {}

for item in taxaSet:
	aValues = []
	if item in hashHits.keys():
		aValues.append(hashHits[item])
	else:
                aValues.append("0")

	if item in hashMisses.keys():
                aValues.append(hashMisses[item])
	else:
                aValues.append("0")
	
	hashOutputTable[item] = aValues


outfile = sys.argv[3] #snakemake.output[0]
handle = open(outfile, "a+") #open handle
handle.write("\t".join(["Taxon", "Fractional_abundance", "Hits_or_misses"]))
handle.write("\n")

for key, array in hashOutputTable.items():
	slicedTaxonomy=";".join(key.split(",")[1:4]).replace("uncultured_marine_bacterium","uncultured").replace("uncultured_bacterium","uncultured")
	if ( float(hashOutputTable[key][0]) >= 0.01 ) or ( float(hashOutputTable[key][1]) >= 0.01 ):
		handle.write("\t".join([slicedTaxonomy, str(hashOutputTable[key][0]), "Matched to Primer"])) 
		handle.write("\n")
		handle.write("\t".join([slicedTaxonomy, str(hashOutputTable[key][1]), "Mismatched to Primer"]))
		handle.write("\n")
