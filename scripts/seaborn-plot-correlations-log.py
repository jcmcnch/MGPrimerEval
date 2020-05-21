#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv
from scipy import stats

sns.set(color_codes=True)

dataset = pd.read_csv(snakemake.input[0], sep="\t")#, names=["cluster", "eASV_frac", "MG_SSU_rRNA_frac"])
axMax = max(list((pd.DataFrame.max(dataset)))[1:]) * 1.05

plot = sns.regplot(x="eASV_frac", y="MG_SSU_rRNA_frac", data=dataset).set_title(snakemake.params[0])
plt.plot([0, 0.8], [0, 0.8], '--', linewidth=2, label="1:1 line")
plt.xlim(0.00001,axMax)
plt.ylim(0.00001,axMax)
plt.legend(loc='upper left', shadow=True)
plt.yscale('log')
plt.xscale('log')

fig = plot.get_figure()

fig.savefig(snakemake.output[0])

#Print information about strength of correlation (can capture as file by redirecting stdout)
#convert to np array, ignoring first column because it is not numeric
nparray = dataset.iloc[:,1:3].to_numpy()
#regress, getting variables out
slope, intercept, r_value, p_value, std_err = stats.linregress(nparray)

#get sample numbers that correspond to SRA ids
aSamples = []
hashSampleSRAConversionTable = {}
for astrLine in csv.reader(open(snakemake.input[1]), csv.excel_tab):
    aSamples.append(astrLine[0])
    hashSampleSRAConversionTable[astrLine[1]] = astrLine[0]
SRAid=str(snakemake.params[0])
BioGEOTRACES_id=str(hashSampleSRAConversionTable[SRAid])

#output handling
with open(snakemake.output[1], "a", newline='\n') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	writer.writerow(["BioGEOTRACES_ID","SRAid","input_file","R^2","p-value","std_err","slope","intercept"])
	writer.writerow([BioGEOTRACES_id,SRAid,str(snakemake.input[0]),str(r_value**2),str(p_value),str(std_err),str(slope),str(intercept)])
