import pandas as pd
import numpy as nu
import os



directory = "/home/lmirete/"
archivo = 'IMMPR11.csv'


#todos = 'Llista completa 323 gens'
#importantes = 'Llista CVID 91 gens'

panel = pd.read_excel('2_Llista gens Panell IDPs NGS VH.xlsx', 'Llista completa 323 gens' , header=None, names = ('SYMBOL','C','C','C','U'))
panel = panel.iloc[:,0].astype(str)


importantes = pd.read_excel('2_Llista gens Panell IDPs NGS VH.xlsx', 'Llista CVID 91 gens' , header=None, names = ('SYMBOL','C','C','C','U'))
importantes = importantes.iloc[:,0].astype(str)

#monogenicos = ['NFKB1','LRBA','CTLA4','PIK3R1','IKZF1','DIKC1','IKBKG']
#TACI = ['TNFRSF13B']
#dudoso = ['BTK?']


diagnostico = pd.read_excel('1_Dades vcf IDCV.xlsx', 'CVID Cohort')
diagnostico = diagnostico[['Panell VH ID','Resultat','Edat','Sexe']]
diagnostico = diagnostico.rename(columns={'Panell VH ID':'SYMBOL','Resultat':'RESULT', 'Edat':'AGE','Sexe':'SEX'})
diagnostico["SEX"] = diagnostico["SEX"].replace({"Masculí": "M", "Femení": "F"})
diagnostico['RESULT'] = diagnostico['RESULT'].replace({'No concloent': 'NC', 'BTK?': 'UNK', 'TNFRSF13B':'TACI','TNFSRF13B':'TACI', 'NFKB1':'MONO','LRBA':'MONO','CTLA4':'MONO','PIK3R1':'MONO','IKZF1':'MONO','DIKC1':'MONO','IKBKG':'MONO', 'DKC1':'MONO'})
diagnostico['DIAGN'] = 'CVID'
diagnostico.set_index('SYMBOL', inplace=True)

fenotipo = pd.read_csv(directory+'fenotipo/fenotips_colobran.csv', delimiter = ' ')
fenotipo = fenotipo[['ID','SEX','AGE_RECRUITMENT']]
fenotipo = fenotipo.rename(columns={'ID':'SYMBOL','AGE_RECRUITMENT':'AGE'})
fenotipo["SEX"] = fenotipo["SEX"].replace({"MALE": "M", "FEMALE": "F"})
fenotipo['DIAGN'] = 'CONTROL'
fenotipo['RESULT'] = 'CONTROL'
fenotipo.set_index('SYMBOL', inplace=True)

diagnostico = pd.concat([diagnostico, fenotipo])



#Variables con ligeras diferencias HH

def df_TODO_HH(archivo, directory, genes):

    #Originar df a partir de la columna proveniente de 'archivo'
    df = pd.read_csv(archivo, delimiter = ',')
    tmp = df[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_count = tmp_groups.count()
    df_panel = pd.merge(tmp_count, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    #Añadir el resto de información al df
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL') 
                tmp_count1 = tmp_groups1.count()
                df_panel1 = pd.merge(df_variations, tmp_count1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})
    
    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0, inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations 

'''
cwd = os.getcwd()  # Get the current working directory (cwd)
files = os.listdir(cwd)  # Get all the files in that directory
print("Files in %r: %s" % (cwd, files))
'''



def df_UNICO_HH(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_unique = tmp_groups.nunique()
    df_panel = pd.merge(tmp_unique, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})


    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL') 
                tmp_unique1 = tmp_groups1.nunique()
                df_panel1 = pd.merge(df_variations, tmp_unique1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0, inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations 



def df_TODO_H(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','freq']]
    df_mask=tmp['freq']=='H'
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_count = tmp_groups.count()
    df_panel = pd.merge(tmp_count, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','freq']]
                df_mask1=tmp1['freq']=='H'
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_count1 = tmp_groups1.count()
                df_panel1 = pd.merge(df_variations, tmp_count1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0, inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations 



def df_UNICO_H(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','freq']]
    df_mask=tmp['freq']=='H'
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_unique = tmp_groups.nunique()
    df_panel = pd.merge(tmp_unique, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','freq']]
                df_mask1=tmp1['freq']=='H'
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_unique1 = tmp_groups1.nunique()
                df_panel1 = pd.merge(df_variations, tmp_unique1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations             



def df_TODO_MISS_HH(archivo, directory,genes):
    
    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','Consequence']]
    df_mask=tmp['Consequence'].str.contains('missense', regex=False)
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_count = tmp_groups.count()
    df_panel = pd.merge(tmp_count, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','Consequence']]
                df_mask1=tmp1['Consequence'].str.contains('missense', regex=False)
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_count1 = tmp_groups1.count()
                df_panel1 = pd.merge(df_variations, tmp_count1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations



def df_UNICO_MISS_HH(archivo, directory,genes):
    
    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','Consequence']]
    df_mask=tmp['Consequence'].str.contains('missense', regex=False)
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_unique = tmp_groups.nunique()
    df_panel = pd.merge(tmp_unique, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','Consequence']]
                df_mask1=tmp1['Consequence'].str.contains('missense', regex=False)
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_unique1 = tmp_groups1.nunique()
                df_panel1 = pd.merge(df_variations, tmp_unique1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations



def df_TODO_MISS_H(archivo, directory,genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','Consequence','freq']]
    df_mask= (tmp['freq']=='H') & (tmp['Consequence'].str.contains('missense', regex=False))
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_count = tmp_groups.count()
    df_panel = pd.merge(tmp_count, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','Consequence','freq']]
                df_mask1= (tmp1['freq']=='H') & (tmp1['Consequence'].str.contains('missense', regex=False))
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_count1 = tmp_groups1.count()
                df_panel1 = pd.merge(df_variations, tmp_count1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations


def df_UNICO_MISS_H(archivo, directory,genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','Consequence','freq']]
    df_mask=(tmp['freq']=='H') & (tmp['Consequence'].str.contains('missense', regex=False))
    tmp = tmp[df_mask]
    tmp = tmp[['SYMBOL','#Uploaded_variation']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_unique = tmp_groups.nunique()
    df_panel = pd.merge(tmp_unique, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','Consequence','freq']]
                df_mask1=(tmp1['freq']=='H') & (tmp1['Consequence'].str.contains('missense', regex=False))
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_unique1 = tmp_groups1.nunique()
                df_panel1 = pd.merge(df_variations, tmp_unique1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations



def df_TODO_POLI_HH(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','PolyPhen']]
    df_mask=tmp['PolyPhen'].str.contains('damaging', regex=False)
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_count = tmp_groups.count()
    df_panel = pd.merge(tmp_count, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','PolyPhen']]
                df_mask1=tmp1['PolyPhen'].str.contains('damaging', regex=False)
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_count1 = tmp_groups1.count()
                df_panel1 = pd.merge(df_variations, tmp_count1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations


def df_UNICO_POLI_HH(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','PolyPhen']]
    df_mask=tmp['PolyPhen'].str.contains('damaging', regex=False)
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_unique = tmp_groups.nunique()
    df_panel = pd.merge(tmp_unique, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','PolyPhen']]
                df_mask1=tmp1['PolyPhen'].str.contains('damaging', regex=False)
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_unique1 = tmp_groups1.nunique()
                df_panel1 = pd.merge(df_variations, tmp_unique1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    # df_variations.loc[:,df_variations.apply(pd.Series.nunique) != 1]
    # para quitar col con varianza 0, no funciona porque yo necesito quitarla en las filas y despues transponer
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations




def df_TODO_POLI_H(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','PolyPhen','freq']]
    df_mask= (tmp['freq']=='H') & (tmp['PolyPhen'].str.contains('damaging', regex=False))
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_count = tmp_groups.count()
    df_panel = pd.merge(tmp_count, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','PolyPhen','freq']]
                df_mask1= (tmp1['freq']=='H') & (tmp1['PolyPhen'].str.contains('damaging', regex=False))
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_count1 = tmp_groups1.count()
                df_panel1 = pd.merge(df_variations, tmp_count1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations



def df_UNICO_POLI_H(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','PolyPhen','freq']]
    df_mask=(tmp['freq']=='H') & (tmp['PolyPhen'].str.contains('damaging', regex=False))
    tmp = tmp[df_mask]
    tmp = tmp[['SYMBOL','#Uploaded_variation']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_unique = tmp_groups.nunique()
    df_panel = pd.merge(tmp_unique, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','PolyPhen','freq']]
                df_mask1=(tmp1['freq']=='H') & (tmp1['PolyPhen'].str.contains('damaging', regex=False))
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_unique1 = tmp_groups1.nunique()
                df_panel1 = pd.merge(df_variations, tmp_unique1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations




def df_TODO_SIFT_HH(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','SIFT']]
    df_mask=tmp['SIFT'].str.contains('deleterious', regex=False)
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_count = tmp_groups.count()
    df_panel = pd.merge(tmp_count, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','SIFT']]
                df_mask1=tmp1['SIFT'].str.contains('deleterious', regex=False)
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_count1 = tmp_groups1.count()
                df_panel1 = pd.merge(df_variations, tmp_count1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations




def df_UNICO_SIFT_HH(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','SIFT']]
    df_mask=tmp['SIFT'].str.contains('deleterious', regex=False)
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_unique = tmp_groups.nunique()
    df_panel = pd.merge(tmp_unique, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','SIFT']]
                df_mask1=tmp1['SIFT'].str.contains('deleterious', regex=False)
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_unique1 = tmp_groups1.nunique()
                df_panel1 = pd.merge(df_variations, tmp_unique1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0,inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations



def df_TODO_SIFT_H(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['SYMBOL','#Uploaded_variation','SIFT','freq']]
    df_mask=(tmp['freq']=='H') & (tmp['SIFT'].str.contains('deleterious', regex=False))
    tmp = tmp[df_mask]
    tmp = tmp[['SYMBOL','#Uploaded_variation']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_count = tmp_groups.count()
    df_panel = pd.merge(tmp_count, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','SIFT','freq']]
                df_mask1=(tmp1['freq']=='H') & (tmp1['SIFT'].str.contains('deleterious', regex=False))
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_count1 = tmp_groups1.count()
                df_panel1 = pd.merge(df_variations, tmp_count1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0, inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations



def df_UNICO_SIFT_H(archivo, directory, genes):
    """
    [summary]

    Args:
        archivo ([.csv]): [.csv, dataframe of one patient]
        directory ([dir]): [directory where are all the patients/controls]
        genes ([Serie]): [d]

    Returns:
        [data frame]: [description]
    """
    

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['SYMBOL','#Uploaded_variation','SIFT','freq']]
    df_mask=(tmp['freq']=='H') & (tmp['SIFT'].str.contains('deleterious', regex=False))
    tmp = tmp[df_mask]
    tmp = tmp[['SYMBOL','#Uploaded_variation']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_unique = tmp_groups.nunique()
    df_panel = pd.merge(tmp_unique, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            # Remove the last 4 letters (.csv)
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','SIFT','freq']]
                df_mask1=(tmp1['freq']=='H') & (tmp1['SIFT'].str.contains('deleterious', regex=False))
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_unique1 = tmp_groups1.nunique()
                df_panel1 = pd.merge(df_variations, tmp_unique1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})
    
    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all', axis=0, inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations 




def df_TODO_REVEL_HH(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','REVEL']]
    df_mask=tmp['REVEL'] >= 0.5
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_count = tmp_groups.count()
    df_panel = pd.merge(tmp_count, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','REVEL']]
                df_mask1=tmp1['REVEL'] >= 0.5
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_count1 = tmp_groups1.count()
                df_panel1 = pd.merge(df_variations, tmp_count1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all',axis=0, inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations




def df_UNICO_REVEL_HH(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['#Uploaded_variation','SYMBOL','REVEL']]
    df_mask=tmp['REVEL'] >= 0.5
    tmp = tmp[df_mask]
    tmp = tmp[['#Uploaded_variation','SYMBOL']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_unique = tmp_groups.nunique()
    df_panel = pd.merge(tmp_unique, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','REVEL']]
                df_mask1=tmp1['REVEL'] >= 0.5
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_unique1 = tmp_groups1.nunique()
                df_panel1 = pd.merge(df_variations, tmp_unique1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all',axis=0, inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations



def df_TODO_REVEL_H(archivo, directory, genes):

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['SYMBOL','#Uploaded_variation','REVEL','freq']]
    df_mask=(tmp['freq']=='H') & (tmp['REVEL'] >= 0.5)
    tmp = tmp[df_mask]
    tmp = tmp[['SYMBOL','#Uploaded_variation']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_count = tmp_groups.count()
    df_panel = pd.merge(tmp_count, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','REVEL','freq']]
                df_mask1=(tmp1['freq']=='H') & (tmp1['REVEL'] >= 0.5)
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_count1 = tmp_groups1.count()
                df_panel1 = pd.merge(df_variations, tmp_count1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})

    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all',axis=0, inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations



def df_UNICO_REVEL_H(archivo, directory, genes):
    """
    [summary]

    Args:
        archivo ([.csv]): [.csv, dataframe of one patient]
        directory ([dir]): [directory where are all the patients/controls]
        genes ([Serie]): [d]

    Returns:
        [data frame]: [description]
    """

    df = pd.read_csv(archivo, delimiter = ",")
    tmp = df[['SYMBOL','#Uploaded_variation','REVEL','freq']]
    df_mask=(tmp['freq']=='H') & (tmp['REVEL'] >= 0.5)
    tmp = tmp[df_mask]
    tmp = tmp[['SYMBOL','#Uploaded_variation']]
    tmp_groups = tmp.groupby('SYMBOL')
    tmp_unique = tmp_groups.nunique()
    df_panel = pd.merge(tmp_unique, genes, on='SYMBOL', how='right')
    df_variations = df_panel.rename(columns={"#Uploaded_variation": str((archivo)[:-4])})

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".csv") and file != archivo:
            name = (str((filename)[:-4])).split('_', 1)[0]
            # Remove the last 4 letters (.csv)
            with open(filename) as readfile:
                df1 = pd.read_csv(readfile, delimiter = ',')
                tmp1 = df1[['#Uploaded_variation','SYMBOL','REVEL','freq']]
                df_mask1=(tmp1['freq']=='H') & (tmp1['REVEL'] >= 0.5)
                tmp1 = tmp1[df_mask1]
                tmp1 = tmp1[['#Uploaded_variation','SYMBOL']]
                tmp_groups1 = tmp1.groupby('SYMBOL')
                tmp_unique1 = tmp_groups1.nunique()
                df_panel1 = pd.merge(df_variations, tmp_unique1, on='SYMBOL', how='left')
                df_variations = df_panel1.rename(columns={"#Uploaded_variation": name})
    
    df_variations.set_index('SYMBOL',inplace=True)
    df_variations.dropna(how='all',axis=0, inplace=True)
    df_variations.dropna(how='all', axis=1, inplace=True)
    df_variations.fillna(0, inplace=True)
    return df_variations     

  

TODO_HH_PANEL = (df_TODO_HH(archivo, directory, panel)).transpose()
TODO_HH_PANEL =TODO_HH_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_HH_PANEL.to_csv(directory+'mis pruebas 3.0/TODO_HH_PANEL.csv')

TODO_HH_IMPORTANTES = (df_TODO_HH(archivo, directory, importantes)).transpose()
TODO_HH_IMPORTANTES =TODO_HH_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_HH_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/TODO_HH_IMPORTANTES.csv')


UNICO_HH_PANEL = (df_UNICO_HH(archivo, directory, panel)).transpose()
UNICO_HH_PANEL =UNICO_HH_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_HH_PANEL.to_csv(directory+'mis pruebas 3.0/UNICO_HH_PANEL.csv')

UNICO_HH_IMPORTANTES = (df_UNICO_HH(archivo, directory, importantes)).transpose()
UNICO_HH_IMPORTANTES =UNICO_HH_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_HH_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/UNICO_HH_IMPORTANTES.csv')

TODO_H_PANEL = (df_TODO_H(archivo, directory, panel)).transpose()
TODO_H_PANEL =TODO_H_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_H_PANEL.to_csv(directory+'mis pruebas 3.0/TODO_H_PANEL.csv')

TODO_H_IMPORTANTES = (df_TODO_H(archivo, directory, importantes)).transpose()
TODO_H_IMPORTANTES =TODO_H_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_H_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/TODO_H_IMPORTANTES.csv')

UNICO_H_PANEL = (df_UNICO_H(archivo, directory, panel)).transpose()
UNICO_H_PANEL =UNICO_H_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_H_PANEL.to_csv(directory+'mis pruebas 3.0/UNICO_H_PANEL.csv')

UNICO_H_IMPORTANTES = (df_UNICO_H(archivo, directory, importantes)).transpose()
UNICO_H_IMPORTANTES =UNICO_H_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_H_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/UNICO_H_IMPORTANTES.csv')


TODO_POLI_HH_PANEL = (df_TODO_POLI_HH(archivo, directory, panel)).transpose()
TODO_POLI_HH_PANEL =TODO_POLI_HH_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_POLI_HH_PANEL.to_csv(directory+'mis pruebas 3.0/TODO_POLI_HH_PANEL.csv')

TODO_POLI_HH_IMPORTANTES = (df_TODO_POLI_HH(archivo, directory, importantes)).transpose()
TODO_POLI_HH_IMPORTANTES =TODO_POLI_HH_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_POLI_HH_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/TODO_POLI_HH_IMPORTANTES.csv')

UNICO_POLI_HH_PANEL = (df_UNICO_POLI_HH(archivo, directory, panel)).transpose()
UNICO_POLI_HH_PANEL =UNICO_POLI_HH_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_POLI_HH_PANEL.to_csv(directory+'mis pruebas 3.0/UNICO_POLI_HH_PANEL.csv')

UNICO_POLI_HH_IMPORTANTES = (df_UNICO_POLI_HH(archivo, directory, importantes)).transpose()
UNICO_POLI_HH_IMPORTANTES =UNICO_POLI_HH_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_POLI_HH_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/UNICO_POLI_HH_IMPORTANTES.csv')

TODO_POLI_H_PANEL = (df_TODO_POLI_H(archivo, directory, panel)).transpose()
TODO_POLI_H_PANEL =TODO_POLI_H_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_POLI_H_PANEL.to_csv(directory+'mis pruebas 3.0/TODO_POLI_H_PANEL.csv')

TODO_POLI_H_IMPORTANTES = (df_TODO_POLI_H(archivo, directory, importantes)).transpose()
TODO_POLI_H_IMPORTANTES =TODO_POLI_H_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_POLI_H_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/TODO_POLI_H_IMPORTANTES.csv')


UNICO_POLI_H_PANEL = (df_UNICO_POLI_H(archivo, directory, panel)).transpose()
UNICO_POLI_H_PANEL = UNICO_POLI_H_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_POLI_H_PANEL.to_csv(directory+'mis pruebas 3.0/UNICO_POLI_H_PANEL.csv')

UNICO_POLI_H_IMPORTANTES = (df_UNICO_POLI_H(archivo, directory, importantes)).transpose()
UNICO_POLI_H_IMPORTANTES =UNICO_POLI_H_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_POLI_H_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/UNICO_POLI_H_IMPORTANTES.csv')


TODO_SIFT_HH_PANEL = (df_TODO_SIFT_HH(archivo, directory, panel)).transpose()
TODO_SIFT_HH_PANEL =TODO_SIFT_HH_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_SIFT_HH_PANEL.to_csv(directory+'mis pruebas 3.0/TODO_SIFT_HH_PANEL.csv')

TODO_SIFT_HH_IMPORTANTES = (df_TODO_SIFT_HH(archivo, directory, importantes)).transpose()
TODO_SIFT_HH_IMPORTANTES =TODO_SIFT_HH_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_SIFT_HH_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/TODO_SIFT_HH_IMPORTANTES.csv')

UNICO_SIFT_HH_PANEL = (df_UNICO_SIFT_HH(archivo, directory, panel)).transpose()
UNICO_SIFT_HH_PANEL =UNICO_SIFT_HH_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_SIFT_HH_PANEL.to_csv(directory+'mis pruebas 3.0/UNICO_SIFT_HH_PANEL.csv')

UNICO_SIFT_HH_IMPORTANTES = (df_UNICO_SIFT_HH(archivo, directory, importantes)).transpose()
UNICO_SIFT_HH_IMPORTANTES =UNICO_SIFT_HH_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_SIFT_HH_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/UNICO_SIFT_HH_IMPORTANTES.csv')

TODO_SIFT_H_PANEL = (df_TODO_SIFT_H(archivo, directory, panel)).transpose()
TODO_SIFT_H_PANEL =TODO_SIFT_H_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_SIFT_H_PANEL.to_csv(directory+'mis pruebas 3.0/TODO_SIFT_H_PANEL.csv')

TODO_SIFT_H_IMPORTANTES = (df_TODO_SIFT_H(archivo, directory, importantes)).transpose()
TODO_SIFT_H_IMPORTANTES =TODO_SIFT_H_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_SIFT_H_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/TODO_SIFT_H_IMPORTANTES.csv')


UNICO_SIFT_H_PANEL = (df_UNICO_SIFT_H(archivo, directory, panel)).transpose()
UNICO_SIFT_H_PANEL = UNICO_SIFT_H_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_SIFT_H_PANEL.to_csv(directory+'mis pruebas 3.0/UNICO_SIFT_H_PANEL.csv')

UNICO_SIFT_H_IMPORTANTES = (df_UNICO_SIFT_H(archivo, directory, importantes)).transpose()
UNICO_SIFT_H_IMPORTANTES = UNICO_SIFT_H_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_SIFT_H_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/UNICO_SIFT_H_IMPORTANTES.csv')


TODO_MISS_HH_PANEL = (df_TODO_MISS_HH(archivo, directory, panel)).transpose()
TODO_MISS_HH_PANEL =TODO_MISS_HH_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_MISS_HH_PANEL.to_csv(directory+'mis pruebas 3.0/TODO_MISS_HH_PANEL.csv')

TODO_MISS_HH_IMPORTANTES = (df_TODO_MISS_HH(archivo, directory, importantes)).transpose()
TODO_MISS_HH_IMPORTANTES =TODO_MISS_HH_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_MISS_HH_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/TODO_MISS_HH_IMPORTANTES.csv')

UNICO_MISS_HH_PANEL = (df_UNICO_MISS_HH(archivo, directory, panel)).transpose()
UNICO_MISS_HH_PANEL =UNICO_MISS_HH_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_MISS_HH_PANEL.to_csv(directory+'mis pruebas 3.0/UNICO_MISS_HH_PANEL.csv')

UNICO_MISS_HH_IMPORTANTES = (df_UNICO_MISS_HH(archivo, directory, importantes)).transpose()
UNICO_MISS_HH_IMPORTANTES =UNICO_MISS_HH_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_MISS_HH_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/UNICO_MISS_HH_IMPORTANTES.csv')

TODO_MISS_H_PANEL = (df_TODO_MISS_H(archivo, directory, panel)).transpose()
TODO_MISS_H_PANEL =TODO_MISS_H_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_MISS_H_PANEL.to_csv(directory+'mis pruebas 3.0/TODO_MISS_H_PANEL.csv')

TODO_MISS_H_IMPORTANTES = (df_TODO_MISS_H(archivo, directory, importantes)).transpose()
TODO_MISS_H_IMPORTANTES =TODO_MISS_H_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_MISS_H_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/TODO_MISS_H_IMPORTANTES.csv')


UNICO_MISS_H_PANEL = (df_UNICO_MISS_H(archivo, directory, panel)).transpose()
UNICO_MISS_H_PANEL = UNICO_MISS_H_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_MISS_H_PANEL.to_csv(directory+'mis pruebas 3.0/UNICO_MISS_H_PANEL.csv')

UNICO_MISS_H_IMPORTANTES = (df_UNICO_MISS_H(archivo, directory, importantes)).transpose()
UNICO_MISS_H_IMPORTANTES = UNICO_MISS_H_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_MISS_H_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/UNICO_MISS_H_IMPORTANTES.csv')


TODO_REVEL_HH_PANEL = (df_TODO_REVEL_HH(archivo, directory, panel)).transpose()
TODO_REVEL_HH_PANEL =TODO_REVEL_HH_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_REVEL_HH_PANEL.to_csv(directory+'mis pruebas 3.0/TODO_REVEL_HH_PANEL.csv')

TODO_REVEL_HH_IMPORTANTES = (df_TODO_REVEL_HH(archivo, directory, importantes)).transpose()
TODO_REVEL_HH_IMPORTANTES =TODO_REVEL_HH_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_REVEL_HH_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/TODO_REVEL_HH_IMPORTANTES.csv')

UNICO_REVEL_HH_PANEL = (df_UNICO_REVEL_HH(archivo, directory, panel)).transpose()
UNICO_REVEL_HH_PANEL =UNICO_REVEL_HH_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_REVEL_HH_PANEL.to_csv(directory+'mis pruebas 3.0/UNICO_REVEL_HH_PANEL.csv')

UNICO_REVEL_HH_IMPORTANTES = (df_UNICO_REVEL_HH(archivo, directory, importantes)).transpose()
UNICO_REVEL_HH_IMPORTANTES =UNICO_REVEL_HH_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_REVEL_HH_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/UNICO_REVEL_HH_IMPORTANTES.csv')

TODO_REVEL_H_PANEL = (df_TODO_REVEL_H(archivo, directory, panel)).transpose()
TODO_REVEL_H_PANEL =TODO_REVEL_H_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_REVEL_H_PANEL.to_csv(directory+'mis pruebas 3.0/TODO_REVEL_H_PANEL.csv')

TODO_REVEL_H_IMPORTANTES = (df_TODO_REVEL_H(archivo, directory, importantes)).transpose()
TODO_REVEL_H_IMPORTANTES =TODO_REVEL_H_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
TODO_REVEL_H_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/TODO_REVEL_H_IMPORTANTES.csv')


UNICO_REVEL_H_PANEL = (df_UNICO_REVEL_H(archivo, directory, panel)).transpose()
UNICO_REVEL_H_PANEL = UNICO_REVEL_H_PANEL.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_REVEL_H_PANEL.to_csv(directory+'mis pruebas 3.0/UNICO_REVEL_H_PANEL.csv')

UNICO_REVEL_H_IMPORTANTES = (df_UNICO_REVEL_H(archivo, directory, importantes)).transpose()
UNICO_REVEL_H_IMPORTANTES = UNICO_REVEL_H_IMPORTANTES.merge(diagnostico, left_index=True, right_index=True, how='inner')
UNICO_REVEL_H_IMPORTANTES.to_csv(directory+'mis pruebas 3.0/UNICO_REVEL_H_IMPORTANTES.csv')
