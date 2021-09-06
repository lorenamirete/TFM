import pandas as pd
import numpy as np
import os

directory = "/home/lmirete/"

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".csv") and filename.startswith('I'):
        name = str(filename)
        with open(filename) as readfile:
            df = pd.read_csv(readfile, delimiter = ',')
            df_mask= (df['Feature_type'] == 'Transcript') & (df['BIOTYPE'] == 'protein_coding') & (df['CANONICAL'] == 'YES') & (df['freq'] != ('0/0'))
            df = df[df_mask]

            df.to_csv(directory+('CSV_FILTRADO/'+name))