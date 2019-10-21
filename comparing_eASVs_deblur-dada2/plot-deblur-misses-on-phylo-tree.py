#!/usr/bin/env python

import csv
import sys

# Import Tree instance and faces module
from ete3 import Tree, faces, TreeStyle

#Get tree from file
t = Tree(sys.argv[1])

def checkValue(number):
    """
    A mini-function for seeing if fractional values are 
    close to parity based on the division of one by the other.
    """
    rounded = round(number, 2)
    if 0.8 <= rounded <= 1.2:
        return "close enough"
    elif rounded < 0.8:
        return "potentially overestimated by deblur"
    elif rounded > 1.2:
        return "potentially underestimated by deblur"

"""
    else:
        if "Kingdom" in splitTax[0]:
             truncTax=splitTax[-1].split(".")[1]
        else:
             for n in range(len(splitTax)):
                 if "uncultured" in splitTax[n]:
                     splitTax.remove(splitTax[n])
             truncTax=splitTax[-1][5:]
"""


#Some code to slice and dice taxonomy labels to make them more readable
hashTax = {}
for astrLine in csv.reader(open(sys.argv[2]), csv.excel_tab):
    splitTax = astrLine[1].split(";")
    if len(splitTax) >= 4:
        if "Kingdom" in splitTax[0]:
             truncTax=splitTax[3].split(".")[1]
        else:
#             for i in range(0, ( len(splitTax) - 2 ) ):
#                 if "uncultured" in splitTax[i]:
#                     splitTax.remove(splitTax[i])
             truncTax=splitTax[3][5:]
    else:
        truncTax=splitTax[-1][5:]

    hashTax[astrLine[0]] = truncTax

#Conversion table to translate between deblur eASV hashes to most similar parent DADA2
hashConversionTable = {}
for astrLine in csv.reader(open(sys.argv[3]), csv.excel_tab):
    try:
        hashConversionTable[astrLine[0]].append(astrLine[1])
    except KeyError:
        hashConversionTable[astrLine[0]] = [astrLine[1]]

#Create a dictionary with information about the effect of biases using DADA2 as a "ground truth" since it's closer to metagenomic reads
#This is sample-specific
hashFracBias = {}
for astrLine in csv.reader(open(sys.argv[4]), csv.excel_tab):
#    if float(astrLine[1].strip()) == 0:
#        for item in hashConversionTable[astrLine[0]]:
#            hashFracBias[item] = "DADA2_undetected"
    if float(astrLine[2].strip()) == 0 and float(astrLine[1].strip()) != 0:
        for item in hashConversionTable[astrLine[0]]:
            hashFracBias[item] = "deblur_undetected"
    elif float(astrLine[2].strip()) != 0 and float(astrLine[1].strip()) != 0:
        for item in hashConversionTable[astrLine[0]]:
            hashFracBias[item] = round( ( float(astrLine[1]) / float(astrLine[2]) ) , 2 ) 

#Now make an array with all the DADA2 ids that will map to the tree
aDADA2hashes = []
l = hashConversionTable.values()
flat_list = [item for sublist in l for item in sublist]
#DADA2 hashes
for key in hashFracBias.keys():
    #if DADA2 hash in sample-specific hits
    if key in flat_list:
            aDADA2hashes.append(key)

#Define tree layout
def mylayout(node):
    if node.is_leaf():
        # text faces support multiline. We add a text face
        # with the whole description of each leaf.
        taxFace = faces.TextFace(hashTax[node.name])
        #taxFace.background.color = "blue"
        # Note that these faces are added in "aligned" mode
        faces.add_face_to_node(taxFace, node, column=0, aligned=True)

        if node.name.strip() in aDADA2hashes:
            descFace = faces.TextFace(hashFracBias[node.name])
            if hashFracBias[node.name] == "deblur_undetected":
                descFace.background.color = "red"
            elif checkValue(hashFracBias[node.name]) == "close enough":
                descFace.background.color = "green"
            elif checkValue(hashFracBias[node.name]) == "potentially overestimated by deblur":
                descFace.background.color = "orange"
            elif checkValue(hashFracBias[node.name]) == "potentially underestimated by deblur":
                descFace.background.color = "blue"
            else:
                descFace.background.color = "white"
            faces.add_face_to_node(descFace, node, column=1, aligned=True)

# And, finally, Visualize the tree using my own layout function

outfile="results/trees/" + sys.argv[4].split("/")[2].split(".")[0] + ".svg"

ts = TreeStyle()
ts.show_leaf_name = False
ts.layout_fn = mylayout
t.render(outfile, w=2500, tree_style=ts)

