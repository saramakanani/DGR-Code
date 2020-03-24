import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt; plt.rcdefaults()

df = pd.read_csv('test_snps.txt', header =None, delim_whitespace = True)
df = df.rename(columns={list(df)[0]:"site_id",list(df)[1]:"mean_freq",list(df)[2]:"mean_depth",list(df)[3]:"site_prev",
                   list(df)[4]:"allele_props",list(df)[5]:"site_type",list(df)[6]:"gene_id",
                   list(df)[7]:"amino_acids",list(df)[8]:"snps"})

## subsets the data 
df2 = df[['site_id','allele_props']]
df_counts = df2.allele_props.str.split("|", expand = True)
df_alleles = df_counts.rename(columns={list(df_counts)[0]:"A_count",list(df_counts)[1]:"C_count",
                                       list(df_counts)[2]:"T_count",list(df_counts)[3]:"G_count"})

## cleaning data: replace something with an empty space
df_alleles['A_count']=df_alleles['A_count'].str.replace('A:','')
df_alleles['C_count']=df_alleles['C_count'].str.replace('C:','')
df_alleles['T_count']=df_alleles['T_count'].str.replace('T:','')
df_alleles['G_count']=df_alleles['G_count'].str.replace('G:','')

cols = df_alleles.columns
df_test = df_alleles.astype(float)
type(np.array(df_test[cols].iloc[1:2,1:2])[0][0])
df_test[cols].eq(float('0.0'))

## output shows the number of values equal to 0 and puts that value into a new column

df_test['equal_to_zero'] = df_test[cols].eq(float('0.0')).sum(axis=1)
df_test['equal_to_zero'].value_counts() 

## generates histogram 
x = ['Monoallelic','Biallelic','Triallelic']
allele_type = df_test['equal_to_zero'].value_counts()

x_pos = [i for i, _ in enumerate(x)]

plt.bar(x_pos, allele_type, color = 'blue')
plt.xlabel("Allele Type")
plt.ylabel("Frequency Across Sites")

plt.xticks(x_pos, x)
plt.show()
