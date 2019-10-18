#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(color_codes=True)

dataset = pd.read_csv(sys.argv[1], sep="\t", names=["deblur_hash", "DADA2_frac", "deblur_frac"])

plot = sns.regplot(x="DADA2_frac", y="deblur_frac", data=dataset).set_title(sys.argv[2])

fig = plot.get_figure()

fig.savefig(sys.argv[3])
