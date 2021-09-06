import pandas as pd
import numpy as np
import os



directory = "/home/lmirete/"


"""

def revel(row):

    # DELETION: chr1 40991136 CT C -> chr1_40991137_T/-
    # DELETION + TWO ALTERNATE ALLELES: chr1 113907204 GTTTTTTTTTTTTTT G,GTTTTTTTT -> chr1_113907205_TTTTTTTTTTTTTT/-/TTTTTTTT
    # INSERTION: chr1 1013466 T TA -> chr1_1013467_-/A
    # INSERTION + TWO ALTERNATE ALLELES: chr1 151817112 C CTGTGTG,CTGTGTGTG -> chr1_151817113_-/TGTGTG/TGTGTGTG
    # INSERTION + DELETION + TWO ALTERNATE ALLELES: chr16 85920239 CT C,CTT -> chr16_85920240_T/-/TT


    if len(row.ref) != len(row.alt):
        #Si la longitud de la fila ref es diferente a la alternativa
        row.grch38_pos += 1
        #entonces le sumamos 1 a la fila posición, así es como está en uploaded variation
        # row.REF[1:] = A lo que hay en REF empezando por la segunda letra (quitamos la posición [0:])
        row.ref = row.ref[1:] if len(row.ref[1:]) > 0 else '-'
        #si es un solo nucleotido nos quedamos con un guión, si tenemos más de uno le quitamos la primera base
        #print(row.ALT)
        row.alt = '/'.join([alt[1:] if len(alt[1:]) > 0 else '-' for alt in row.alt.split(',')])
        #la función join quiere decir 'unir con lo que está antes del punto a la izquierda'
        #print(row.ALT)
        #actuamos con la primera base igual que actuamos con REF, además eliminamos las posibles comas del nombre

    return 'chr' + str(row['#chr']) + '_' + str(row.grch38_pos) + '_' + str(row.ref) + '/' + str(row.alt)
"""

#revel_scores = "new_tabbed_revel_grch38.tsv"
revel_scores = 'Revel_variation.csv'

def añadir_revel(directory):

#abrir csv records REVEL
    df = pd.read_csv(revel_scores, delimiter = ',', dtype={'#chr': str})
    #df['#Uploaded_variation'] = df.apply(lambda row: revel(row), axis=1)
            #ESTA ÚLTIMA LÍNEA ES MUY ÚTIL E IMPORTANTE, ESTOY ITERANDO EN UNA COL O FILA ESPECÍFICA
            #apply -> para evitar un for, axis -> 0 or ‘index’: apply function to each column. 1 or ‘columns’: apply function to each row.
            #lamba es una funcion anonima. El contenido de una función lambda debe ser una única expresión en lugar de un bloque de acciones

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if (filename.endswith(".csv")) & (filename.startswith('I')) :            
            with open(filename) as readfile:
                csv_df = pd.read_csv(readfile, delimiter = ',')
                #csv_df = df.merge(df[['#Uploaded_variation', 'REVEL']], how='left', on = '#Uploaded_variation')
                df_final = pd.merge(csv_df, df[['#Uploaded_variation', 'REVEL_x']], on='#Uploaded_variation', how='left')

            df_final.to_csv(directory+'INMUNO CSV CON REVEL/'+str(filename), index = False)

            
            #Ver a qué nivel tengo que poner el "guardado"

            # TEST check all #Uploaded_variation have a variation id
            #if not txt_df[txt_df.freq.isnull()].empty:
            #    missing_variationid = txt_df[txt_df.variation_id.isnull()]['#Uploaded_variation'].unique()
            #    print('ERROR: Uploaded_variation missing: ' + ', '.join(missing_variationid))
                
            #txt_df.to_csv(name+'.csv', index=False)
                    
df_cosa = añadir_revel(directory)
#print(df_cosa)
#df_cosa.to_csv(str(filename), index = False)