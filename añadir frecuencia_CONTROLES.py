import pandas as pd
import os


directory = "/home/lmirete/"
archivo = 'IDP0019.txt'


def añadir_frecuencia(directory):
    for file in os.listdir(directory):
        filename = os.fsdecode(file)

        if filename.endswith(".vcf"):
            name = str((filename)[:-4]).split('_', 1)[0]
            # print(name)
            with open(filename) as readfile:
                vcf_df = pd.read_csv(readfile, delimiter = '\t')
                vcf_df['freq'] = vcf_df[(str(name))]
                vcf_df['freq'] = vcf_df['freq'].replace({'1/1': 'H', '0/1': 'HH', '1/2' : 'Hh'})
                vcf_df['#Uploaded_variation'] = vcf_df['ID']
                
            filename = filename.replace('.vcf', '.txt')
            with open(filename) as readfile:
                txt_df = pd.read_csv(readfile, delimiter = '\t')
                txt_df = txt_df.merge(vcf_df[['#Uploaded_variation', 'freq']], how='left', on = '#Uploaded_variation')

            #   TEST check all #Uploaded_variation have a variation id
            if not txt_df[txt_df.freq.isnull()].empty:
                missing_variationid = txt_df[txt_df.variation_id.isnull()]['#Uploaded_variation'].unique()
                print('ERROR: Uploaded_variation missing: ' + ', '.join(missing_variationid))
                
            txt_df.to_csv(name+'.csv', index=False)
                    
añadir_frecuencia(directory)
