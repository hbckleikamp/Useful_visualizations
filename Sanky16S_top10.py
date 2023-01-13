#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 16:03:55 2021

@author: hugokleikamp
"""
#%% clear variables and console

try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

#%% Modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import random, re, requests
import threading, time, string
import pickle


from Bio.SeqUtils import molecular_weight
from itertools import chain, groupby
from collections import Counter
from openpyxl import load_workbook 

import urllib, zipfile

import math, statistics
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


import numpy as np
from scipy.special import comb

#import seaborn as sns; sns.set()
# from sklearn.datasets import load_digits
# from sklearn.decomposition import PCA


#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
print(os.getcwd())



def repunk(x):
    x=x.strip()
    if "unknown" in x:
        x=""
    return x


def Ca_rem(x): #remove candidatus
    if type(x)==str:
        x=x.replace("Ca_","")
        x=x.replace("Candidatus_","")
        
    return x

def smoothstep(x, x_min=0, x_max=1, N=1):
    x = np.clip((x - x_min) / (x_max - x_min), 0, 1)

    result = 0
    for n in range(0, N + 1):
         result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n

    result *= x ** (N + 1)

    return result


class rectangle:
    def __init__(self,taxa,count,color):
        self.taxa=taxa
        self.count=count
        self.color=color #good/dump/gap/missing
        
        #add later:
            # -color (based on taxa and status)
            # -xcor
            # -top,bottom,centre
            # -links




#%% load zotus
headers=dict(pd.read_excel(
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/Data_resources/Experimental_data/16S/Novogene info/DNAsampleheaders.xlsx"
                           ).values.tolist())

samples=list(headers.values())


zotus=pd.read_csv("/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/1_Annotation/16S/1_ASV_picking/zotus.tsv",sep="\t")
zotus=zotus[['#OTU ID']+list(headers.keys())]
zotus.columns=(['Feature ID']+list(headers.values()))


#zotus=zotus.applymap(lambda x: (denoise_z(x)))
zotus["Feature ID"]=zotus["Feature ID"].apply(lambda x: int(x[4:]))   
#for now select a single collumn to count by
countcols=["DXP 2,0","GW R1 2,0","SP 2,0"]
countcol="Count"
zotus=zotus[["Feature ID"]+countcols]
zotus["Count"]=zotus[countcols].sum(axis=1)
zorus=zotus[["Feature ID","Count"]]


#%%


dbs=["GTDB16S_all_FL_hom_tax","MiDAS","Silva"]
    
Dump_keywords=["uncultured","unidentified","organism","Unknown","_metagenome","unknown"  
               "Subgroup", "_group","_bacterium", "_proteobacterium",
               "s__anaerobic_digester","s__human_gut" #custom blacklist 
               
               #"_clade",  
               #"_soil", "_sediment", "_marine","rock_porewater", "_gut", "_forest" 
               #"anaerobic_digester", "toluene-degrading_methanogenic"
               ]              

status=["good","dump","gap","missing"]




files=["/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/Annotation/16S/3_Training_taxonomy/output/GTDB_all_all.tsv",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/Annotation/16S/3_Training_taxonomy/output/GTDB_all_FL_taxonomy.tsv",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/Annotation/16S/3_Training_taxonomy/output/GTDB_reps_all_taxonomy.tsv",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/Annotation/16S/3_Training_taxonomy/output/GTDB_reps_FL_taxonomy.tsv",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/Annotation/16S/3_Training_taxonomy/output/Midas_taxonomy.tsv",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/Annotation/16S/3_Training_taxonomy/output/Silva_taxonomy.tsv"]


files=[
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/1_Annotation/16S/3_Training_taxonomy/output/GTDB_reps_FL_taxonomy.tsv",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/1_Annotation/16S/3_Training_taxonomy/output/Midas_taxonomy.tsv",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/1_Annotation/16S/3_Training_taxonomy/output/Silva_taxonomy.tsv"]

#tranks=["phylum_name","genus_name"]
ranks=["superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"] 
tranks=ranks[1:]

#%% load homonyms

hfiles=["/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/Taxonomy mapping/ssu_silva_taxonomy__to_GTDB_mapping.xlsx",
        "/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/Taxonomy mapping/ncbi_taxonomy__to_GTDB_mapping.xlsx"]

#16S mapping: Silva only
hfiles=[
        "/Volumes/Seagate_SSD/paper 3/third_paper_newest/1_Database_comparison/Taxonomy mapping/ssu_silva_taxonomy__to_GTDB_mapping.xlsx"]

hdfs=[]
for h in hfiles:
    hdf=pd.read_excel(h)
    hdf=hdf.iloc[:,[1,2,-1]]

    hdf=hdf[hdf['gtdb_taxonomy_taxfreq']>0.75]
    hdf.columns=["other_db","gtdb","gtdb_taxonomy_taxfreq"]
    hdf=hdf.set_index("other_db")
    

    hdfs.append(hdf["gtdb"])

hdfs=pd.concat(hdfs)
#%%

  
taxcols=[]
barcols=[]
once=True
for file in files:

    xlsdf=pd.read_csv(file,sep="\t").iloc[1:].reset_index(drop=True).fillna("")
    xlsdf.pop("Confidence")         
    xlsdf["Taxon"]=xlsdf["Taxon"].apply(lambda x: x+";"*(len(ranks)-x.count(";")-1))
    xlsdf[ranks]= xlsdf["Taxon"].str.split(";", expand=True).iloc[:,:len(ranks)]  
    xlsdf=xlsdf.applymap(lambda x: repunk(x))            
    

    #replace homonyms v2
    for target in tranks:

        q=xlsdf[target].isin(hdfs.index)        
        xlsdf.loc[q,target]=hdfs.loc[xlsdf.loc[q,target].values].values
    
    

    xlsdf=xlsdf[['Feature ID']+ ranks]
    bartax=["_".join(Path(file).stem.split("_")[:-1])+'_'+r for r in ranks]
    barcols.append(bartax)
    taxcols.extend(bartax)
    xlsdf.columns=['Feature ID']+ bartax
    
    

    if once:
        mdf=xlsdf
        once=False
    else:
        mdf=mdf.merge(xlsdf,how="left",on="Feature ID")


mdf["Feature ID"]=mdf["Feature ID"].apply(lambda x: int(x[4:]))    
barcols=barcols[::-1]

#add zotu counts
mdf=mdf.merge(zotus[["Feature ID",countcol]],on="Feature ID")


#only for main figure (groupby first then regex)
for i in ["_A","_B","_C","_D","_E","_F","_G","_H"]: #GTDB subspecifications
    mdf[taxcols]=mdf[taxcols].replace(r"\{0}$".format(i), "",regex=True)

#retain in both version
mdf[taxcols]=mdf[taxcols].applymap(lambda x: x.replace(" ","_"))
mdf[taxcols]=mdf[taxcols].applymap(lambda x: x.replace("Ca_",""))
mdf[taxcols]=mdf[taxcols].applymap(lambda x: x.replace("Candidatus_",""))

    
mdf=mdf.groupby(taxcols)["Count"].sum().reset_index()
mdf[countcol]=mdf[countcol]/mdf[countcol].sum()*100 #normalize counts
    
    #color based on phylum
    
    # phyla_cols=[t for t in mdf.columns if "phylum" in t]
    # unitax=np.unique(mdf[phyla_cols].values)

#%% print top 10 taxa

cs=[]    
for i in mdf.columns[:-1]:
    c=mdf[[i,"Count"]].groupby(i)["Count"].sum()
    c=c[c.index!=""] #remove unannotated
                     #remove dump
                     
                     
    dump=pd.concat([c[c.index.str.contains(d)] for d in Dump_keywords])
    dump=dump[~dump.index.duplicated(keep='first')]
    c=c[~c.index.isin(dump.index)].reset_index()
    c.columns=[i,"Count_"+i]
    c=c[~c[i].str.contains("_gap_")]
    cs.append(c)                

cs=pd.concat(cs,axis=1)
cs.to_csv("top1016S.csv")