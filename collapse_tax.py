#!/usr/bin/env python

import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep='\t', index_col=0, header=0)

df = df.fillna('-').max(axis=1, numeric_only=False)
print (df.shape)
df.index.name = 'feature'
df.to_csv('taxonomy.txt', sep = '\t', header = ['lineage'])

