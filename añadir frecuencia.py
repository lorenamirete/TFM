import pandas as pd
import os


directory = "/home/lmirete/"
archivo = 'IMMPR11.txt'


def variant_vcf2txt(row):

    # DELETION: chr1 40991136 CT C -> chr1_40991137_T/-
    # DELETION + TWO ALTERNATE ALLELES: chr1 113907204 GTTTTTTTTTTTTTT G,GTTTTTTTT -> chr1_113907205_TTTTTTTTTTTTTT/-/TTTTTTTT
    # INSERTION: chr1 1013466 T TA -> chr1_1013467_-/A
    # INSERTION + TWO ALTERNATE ALLELES: chr1 151817112 C CTGTGTG,CTGTGTGTG -> chr1_151817113_-/TGTGTG/TGTGTGTG
    # INSERTION + DELETION + TWO ALTERNATE ALLELES: chr16 85920239 CT C,CTT -> chr16_85920240_T/-/TT


    if len(row.REF) != len(row.ALT):
        #Si la longitud de la fila ref es diferente a la alternativa
        row.POS += 1
        #entonces le sumamos 1 a la fila posición, así es como está en uploaded variation
        # row.REF[1:] = A lo que hay en REF empezando por la segunda letra (quitamos la posición [0:])
        row.REF = row.REF[1:] if len(row.REF[1:]) > 0 else '-'
        #si es un solo nucleotido nos quedamos con un guión, si tenemos más de uno le quitamos la primera base
        row.ALT = '/'.join([alt[1:] if len(alt[1:]) > 0 else '-' for alt in row.ALT.split(',')])
        #con join unimos con lo que está antes del punto a la izquierda
        #actuamos con la primera base igual que actuamos con REF, además eliminamos las posibles comas del nombre

        
    return row['#CHROM'] + '_' + str(row.POS) + '_' + row.REF + '/' + row.ALT


def añadir_frecuencia(directory):
    for file in os.listdir(directory):
        filename = os.fsdecode(file)

        if filename.endswith(".vcf"):
            name = str((filename)[:-4]).split('_', 1)[0]
            with open(filename) as readfile:
                vcf_df = pd.read_csv(readfile, delimiter = '\t')
                vcf_df['freq'] = vcf_df.iloc[:,-1].str.split(',').str.get(0).str.split(':').str.get(0)
                vcf_df['freq'] = vcf_df['freq'].replace({'1/1': 'H', '0/1': 'HH', '1/2' : 'Hh'})
                vcf_df['#Uploaded_variation'] = vcf_df.apply(lambda row: variant_vcf2txt(row), axis=1)
                #ÚLTIMA LÍNEA ITERANDO EN UNA COL O FILA ESPECÍFICA
                #apply axis-> 0 or ‘index’: apply function to each column. 1 or ‘columns’: apply function to each row.
                # lamba es una funcion anonima. El contenido de una función lambda debe ser una única expresión en lugar de un bloque de acciones

            filename = filename.replace('.vcf', '.txt')
            with open(filename) as readfile:
                txt_df = pd.read_csv(readfile, delimiter = '\t')
                txt_df = txt_df.merge(vcf_df[['#Uploaded_variation', 'freq']], how='left', on = '#Uploaded_variation')

            # TEST check all #Uploaded_variation have a variation id
            if not txt_df[txt_df.freq.isnull()].empty:
                missing_variationid = txt_df[txt_df.variation_id.isnull()]['#Uploaded_variation'].unique()
                print('ERROR: Uploaded_variation missing: ' + ', '.join(missing_variationid))
                
            txt_df.to_csv(name+'.csv', index=False)
                    
    df_frecuencia = añadir_frecuencia(directory)
    df_frecuencia.to_csv((str((filename)[:-4])).split('_', 1)[0]+'.csv', index = False)

añadir_frecuencia(directory)
