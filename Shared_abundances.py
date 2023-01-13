#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 18:03:57 2021

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

import math, statistics
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns; sns.set()
from sklearn.datasets import load_digits
from sklearn.decomposition import PCA


#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
print(os.getcwd())


#%%
#import venn2,venn3
folder="/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/main_ouput/"
ranks=["superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]  
plants=["DXP","GW","SP"]
methods=["MP","MG","16S"]

from matplotlib_venn import venn3,venn2
import itertools

# VENN diagram
#test comparison with GW DXP SP, GTDB
#ranks=["phylum_name","class_name","order_name","family_name","genus_name","species_name"]  
#ranks=["genus_name"]
labels=["DXP","GW","SP"]

files=os.listdir(folder)
files=[file for file in files if file[0].isalnum() and file.endswith(".xlsx")]


alldf=pd.DataFrame()
for plant in plants:
    
    plantfiles=[file for file in files if plant in file]
    
    for rank in ranks:
        sets=[]
        for file in plantfiles:
            
            
        
            filepath=str(Path(folder,file))
            xlsdf=pd.read_excel(filepath)
            xlsdf=xlsdf[xlsdf[rank+"_count"]>0]
            sets.append(set(xlsdf[rank].dropna().tolist()))
        
        
        fig = plt.figure()
        venn3(sets,methods)
        plt.title(plant+" "+rank.split("_")[0])
        
        
        #get shared %
        aMP_MG_16S = set.intersection(*sets)
        
        aMP_MG =  [i for i in set.intersection( sets[0],sets[1]) if i not in sets[2]]
        aMP_16S = [i for i in set.intersection( sets[0],sets[2]) if i not in sets[1]]
        aMG_16S = [i for i in set.intersection( sets[1],sets[2]) if i not in sets[0]]        
        

    
        contributions=[]
        for file in plantfiles:
            filepath=str(Path(folder,file))
            xlsdf=pd.read_excel(filepath).fillna(0)
            
            for i in [aMP_MG_16S,aMP_MG,aMP_16S,aMG_16S]:
                
                
                ratio=sum(xlsdf.loc[xlsdf[rank].isin(i), rank+"_count"])/sum(xlsdf[rank+"_count"])
                contributions.append(ratio)
                
    
    
        arr=np.array(contributions).reshape(-1,4)
        df=pd.DataFrame(arr,index=methods,columns=["MP_MG_16S","MP_MG","MP_16S","MG_16S"])
        df["self"]=1-df.sum(axis=1)
        df.index=[plant+" "+rank.split("_")[0]+" "+i for i in df.index]
        
            
        #stacked bar
        # fig = plt.figure()
        # df.plot(kind='bar', stacked=True, figsize=(10, 6), legend=False)
        
        colors=np.array([["#d2d2d3","#e7ddcd","#e6cee1","#c6e7e7","#e7b8b9"],
                         ["#d2d2d3","#e7ddcd","#e6cee1","#c6e7e7","#bbdcc2"],
                         ["#d2d2d3","#e7ddcd","#e6cee1","#c6e7e7","#b8c6df"]])
        #sns.set_style("white")
        arr=df.values
        fig, ax = plt.subplots()
        
        xcors=np.linspace(0,0.25*3,3)
        
        for x,i in enumerate(arr):
            bottom=0
            for y,j in enumerate(i):
                
                
        
                ax.bar(xcors[x], height=j, bottom=bottom, width=0.35, color=colors[x,y])
                bottom+=j
                
        plt.axis('off')
        plt.grid(b=None)
                

        plt.show()
                
    
        alldf=alldf.append(df)
            
#%%
alldf.to_excel("shared_abundances.xlsx")
