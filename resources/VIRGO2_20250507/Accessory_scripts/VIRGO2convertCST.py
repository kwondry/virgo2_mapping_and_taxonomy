#!/usr/bin/python
import pandas as pd 
import sys
import numpy as np

#Reading in dataframe, expects first column be 'UID', and renaming columns to be taxa with values of either relative abundances or read counts
data = pd.read_csv(sys.argv[1],sep=',',index_col=0)
data=data.T
data=data.reset_index()
data=data.rename(columns={'index':'sampleID'})
print(data)
#These taxa are being collated to the genus level to conform with expectations of Valencia, which has some CSTs that are typified by a genus and not a species
taxa2b_collated = ['Gardnerella','Fannyhessea']
new_name_dict = {'Gardnerella':'g_Gardnerella','Bifidobacterium':'g_Bifidobacterium','Streptococcus':'g_Streptococcus',
                 'Enterococcus':'g_Enterococcus','Staphylococcus':'g_Staphylococcus','Fannyhessea':'g_Fannyhessea'}

#collating species within target genera
for taxa in taxa2b_collated:
    target_cols = [col for col in data.columns if taxa in col]
    revised_name = new_name_dict[taxa]
    data[revised_name]=data.loc[:,target_cols].sum(axis=1)
    data = data.drop(target_cols,axis=1)

#manually fixing a couple taxa that were split into multiple species since the publication of Valencia 
#fixing prevotella timonensis, which was 
Ptim_cols = ['Prevotella_timonensis','Prevotella_spNov3','Prevotella_spNov2','Prevotella_sp000479005']
data['Prevotella_timonensis']=data.loc[:,Ptim_cols].sum(axis=1)
data = data.drop(Ptim_cols[1:],axis=1)

#fixing speciated Lactos
data['Lactobacillus_gasseri'] = data['Lactobacillus_gasseri'] + data['Lactobacillus_paragasseri']
data = data.drop(['Lactobacillus_paragasseri'],axis=1)
data=data.rename({'UBA629_sp005465875': 'Ca_Lachnocurva_vaginae', 'KA00274_sp902373515': 'Amygdalobacter_indicium','g_Gardnerella':'Gardnerella_vaginalis','g_Fannyhessea':'Fannyhessea_vaginae'}, axis=1)

data=data.set_index('sampleID')

fileOUT = sys.argv[1].split('.')[0]

#sorting dataframe by column sum (taxa sorted by study wide relative abundance)
s = data.sum()
data = data[s.sort_values(ascending=False).index]
data.to_csv("%s_revised4CST.csv" %(fileOUT))