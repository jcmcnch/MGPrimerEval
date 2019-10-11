#!/usr/bin/env python3

import csv, re
from subprocess import check_output
import networkx as nx
from networkx.algorithms.components.connected import connected_components as cc

hashSampleSRAConversionTable = {}

#Conversion table between SRAid and sample ID
for astrLine in csv.reader(open(snakemake.input[0]), csv.excel_tab):
	hashSampleSRAConversionTable[astrLine[1]] = astrLine[0]

hashASVTable = {}

#Convert ASV table to hash
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

hashMG = {}

#Populate dictionary from blastn results
for astrLine in csv.reader(open(snakemake.input[2]), csv.excel_tab):
	try:
		hashMG[astrLine[0]].append(astrLine[1])
	except KeyError:
		hashMG[astrLine[0]] = [astrLine[1]]

#Use networkx to convert dictionary into graph object
G = nx.from_dict_of_lists(hashMG)

hashMembership = {}
#Use connected_components to connect SRA IDs with ASV ids

"""
Note that what we're doing here is comparing short-read SSU rRNA from MG with ASVs.
Since the same MG read could perfectly match multiple ASVs, we're merging them into a cluster. 
Since this cluster will probably not be stable across samples, it is calculated de novo for each sample.
The abundance of ASVs for this cluster are then summed and compared to the fractional abundance of
MG SSU rRNA reads that match this cluster. Clusters are NOT related to sequence identity (i.e. an OTU)
but simply how many overlapping blast hits there are per cluster. We lose a lot of resolution with this
method, but it does give us a rough approximation of the quantitativeness of both methods, and if
there are any gross biases.
"""
for item in list(cc(G)):

	aSRR, aASV = [], [] #Clear lists

	for element in item:

		#if string starts with SRR, add to separate list
		if re.search('^SRR', element):
			aSRR.append(element)

		#else, it's an ASV id so add it to a separate list
		else:
			aASV.append(element)

	#Populate hash with tuple as key
	hashMembership[tuple(aASV)] = aSRR

hashCounts = {}

SRAid=snakemake.params[0]
fileLocation='compareTags2MG/01-concatenated/{}.PROK.concat.515Y-926R.fasta'.format(SRAid)
totalSeqs=check_output(['grep', '-c', 'SRR', fileLocation])
totalSeqs=int(totalSeqs.decode('UTF-8').strip())

sampleID=hashSampleSRAConversionTable[SRAid]

outfile = snakemake.output[0]
handle = open(outfile, "a+") #open handle

#Populate hashCounts with proportions from MG and ASV table
for key, array in hashMembership.items():
	iFracSumASV=0 #the sum of all fractional abundances
	for asvHash in key:
		iLocation = hashASVTable["#OTU ID"].index(asvHash) #Get location in list
		iFracSumASV += float(hashASVTable[sampleID][iLocation]) #Use location to get fraction, add
	iFracSumMG = len(array) / float(totalSeqs) #Divide number of MG hits by total to get fraction
	handle.write("\t".join([";".join(key), str(iFracSumASV), str(iFracSumMG)]))
	handle.write("\n") 
