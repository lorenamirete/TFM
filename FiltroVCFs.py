import tabix
import os
import pandas as pd
import gzip

directorio = '/home/lmirete'
posiciones = pd.read_csv('/home/lmirete/Downloads/posiciones.csv', delimiter = ',')


for file in os.listdir(directorio):
    filename = os.fsdecode(file)
    if filename.endswith(".vcf.gz") and filename.startswith("I"):
        
        # Extracción vcf records
        records = []
        tb=tabix.open(filename)
        for index, row in posiciones.iterrows():
            records.extend(tb.query(row.chromosome_name, row.start_position, row.end_position))

        # Extracción header
        with gzip.open(filename) as fh:
            header = fh.read().decode('utf-8').split('\n')[:68]

        # Escribir header + records en formato vcf
        name = str(filename)[:-7]
        with open(name + '_filtered.vcf', 'w') as fw:
            fw.write('\n'.join(header) + '\n')
            for record in records:
                fw.write('\t'.join(record) + '\n')
