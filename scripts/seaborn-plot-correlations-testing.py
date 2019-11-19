#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

sns.set(color_codes=True)

dataset = pd.read_csv(sys.argv[1], sep="\t", names=["cluster", "eASV_frac", "MG_SSU_rRNA_frac"])
axMax = max(list((pd.DataFrame.max(dataset)))[1:]) * 1.05
plot = sns.regplot(x="eASV_frac", y="MG_SSU_rRNA_frac", data=dataset).set_title(str(sys.argv[1]))
plt.plot([0, 0.8], [0, 0.8], '--', linewidth=2, label="1:1 line")
plt.xlim(-0.02,axMax)
plt.ylim(-0.02,axMax)
plt.legend(loc='upper left', shadow=True)
#(xlim=(0,10),ylim=(0,100)))
#set(xlim=(0,15),ylim=(0,100))

fig = plot.get_figure()

fig.savefig(str(sys.argv[1]) + ".svg")
