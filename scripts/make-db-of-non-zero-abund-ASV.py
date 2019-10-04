#!/usr/bin/env python3

import sys
import argparse
import csv

aSamples = []
hashSampleSRAConversionTable = {}

#Use conversion table
for astrLine in csv.reader(open(snakemake.input[0]), csv.excel_tab):
    aSamples.append(astrLine[0])
    hashSampleSRAConversionTable[astrLine[1]] = astrLine[0]

hashASVTable = {}

with open(snakemake.input[1], newline="") as f:
    reader = csv.reader(f, csv.excel_tab)
    aHeaders = next(reader) #Get headers from first line

    #populate dictionary keys
    for item in aHeaders:
        hashASVTable[item.strip()] = []

    #populate array with rest
    for row in reader:
        for i in range(0,len(aHeaders)):
            hashASVTable[aHeaders[i]].append(row[i].strip())

key = hashSampleSRAConversionTable[snakemake.params[0]] #get sample ID we want to analyze
aAbundances = hashASVTable[key]

for i in range(0, len(aAbundances)): #iterate across rows
    outfile = snakemake.output[0]
    handle = open(outfile, "a+") #open handle
    if float(aAbundances[i]) > 0: #only record ASV id if abundance > 0
        handle.write(hashASVTable["#OTU ID"][i] + "\n") #append to file
