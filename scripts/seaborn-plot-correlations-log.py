#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(color_codes=True)

dataset = pd.read_csv(snakemake.input[0], sep="\t")

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
