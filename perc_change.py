#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 10:53:17 2021

@author: hugokleikamp
"""
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

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
sns.set_style("white")

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
print(os.getcwd())


#%% fold change
ranks=["superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]  
plants=["DXP","GW","SP"]
methods=["MP","S16","MG"]

class dataset:
        def __init__(self,name,rank,counts,method,plant):
            self.name=name
            self.rank=rank
            self.counts=counts
            self.method=method
            self.plant=plant
            
           
      
folders=["/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/main_ouput/"]    
for folder in folders:
            
    #folder="/Volumes/Seagate_SSD/third_paper new/Annotation/Output files/rank_quant/"
    
    
    datasets=[]             
    for file in os.listdir(folder):
        if file[0].isalnum() and file.endswith(".xlsx"):
            
#
            
            plant=[i for i in plants if i in file][0]
            method=[i for i in methods if i in file][0]
            
    
            
            xlsdf=pd.read_excel(str(Path(folder,file)))
            
            for rank in ranks:
                name=" ".join([rank.split("_")[0],plant,method])
                
                counts=xlsdf[[rank,rank+"_count"]]
                counts=counts[~counts[rank].isnull()]
                counts=counts.rename(columns={name:(plant+" "+method)})
                mydata=dataset(name,rank,counts,method,plant)
                datasets.append(mydata)
    
    
    
    
    #%% get Biases
    
    
    once=True
    
    for plant in plants:    
        
    
        
    
        allfolds=pd.DataFrame()
        
        for rank in ranks:
        #%%
            first1=True
            names=[]
            combdf=pd.DataFrame()
            
            for data in datasets:
        #%%
                if data.rank==rank and data.plant==plant:
                    
                    
        
                    names.append(data.name)
                    
                    if first1:
                        first1=False
                        combdf=data.counts
                    else:
                        combdf=combdf.merge(data.counts,on=data.rank,how='outer') #how=left
                        
            combdf=combdf.set_index(rank,drop=True).fillna(0)
            combdf.columns=[" ".join(name.split(" ")[1:]) for name in names]
            
        
            MGcol=[i for i in combdf.columns if "MG" in i][0]
            MPcol=[i for i in combdf.columns if "MP" in i][0]
            S16col=[i for i in combdf.columns if "S16" in i][0]
            
            foldnames=["MG-MP perc_change",
                       "MP-MG perc_change",
                       "16S-MP perc_change", 
                       "MP-16S perc_change", 
                       "16S-MG perc_change",
                       "MG-16S perc_change"]
            
            
            folds=[[MPcol,MGcol], 
                   [MGcol,MPcol], 
                   [MPcol,S16col], 
                   [S16col,MPcol], 
                   [MGcol,S16col],
                   [S16col,MGcol]]
            
            for ix,f in enumerate(folds):
           
                comp=combdf[f].apply(lambda x: (x.iloc[0]-x.iloc[1])/((x.iloc[0]+x.iloc[1])/2),axis=1)
                comp.name=foldnames[ix]+" "+plant#+" "+change
                combdf=combdf.merge(comp,how="left",on=rank)
            
            allfolds=allfolds.append(combdf)
            #%%
        allfolds.index.name = 'taxa'
        
    
        if once==True:
            allplants=allfolds
            once=False
        else:
            allplants=allplants.merge(allfolds,on="taxa",how='outer')
    
    
    
    
    allplants.to_excel("perc_change.xlsx")