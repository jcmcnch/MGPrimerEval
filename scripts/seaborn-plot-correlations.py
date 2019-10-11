#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set(color_codes=True)

dataset = pd.read_csv(snakemake.input[0], sep="\t", names=["cluster", "eASV_frac", "MG_SSU_rRNA_frac"])

plot = sns.regplot(x="eASV_frac", y="MG_SSU_rRNA_frac", data=dataset).set_title(snakemake.params[0])

fig = plot.get_figure()

fig.savefig(snakemake.output[0])
