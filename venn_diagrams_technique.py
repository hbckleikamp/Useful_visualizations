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
import seaborn as sns; sns.set()
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
#folder="/Volumes/Seagate_SSD/third_paper new/Annotation/Output files/rank_quant/"

folders=["/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/main_ouput/"]

#folders=["/Volumes/Seagate_SSD/third_paper_newest/1_Database_comparison/2_Quantification/outputs/blca/"]

for folder in folders:

    
    ranks=["superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]  
    plants=["DXP","GW","SP"]
    methods=["MP_quant","MG_quant","S16_quant"]
    
    from matplotlib_venn import venn3,venn2
    import itertools
    
    # VENN diagram
    #test comparison with GW DXP SP, GTDB
    ranks=["phylum_name","class_name","order_name","family_name","genus_name","species_name"]  
    #ranks=["genus_name"]
    labels=["DXP","GW","SP"]
    
    files=os.listdir(folder)
    files=[file for file in files if file[0].isalnum() and file.endswith(".xlsx")]
    
    
    
    mlabels=["MP","MG","16S"]
    
    
        #for plant in plants:
            
    
    for ix,plant in enumerate(plants):
        
        alls=[]
        means=[]
        stds=[]
        for rank in ["genus_name"]:
        
            plantfiles=[file for file in files if plant in file and file[0].isalnum()]
            
            
            methodfiles=[p for m in methods for p in plantfiles  if m in p and p[0].isalnum()]

        
            sets=[]
            for file in methodfiles:
                
                
            
                filepath=str(Path(folder,file))
                xlsdf=pd.read_excel(filepath)
                xlsdf=xlsdf[xlsdf[rank+"_count"]>0]
                sets.append(set(xlsdf[rank].dropna().tolist()))
            
        
            fig = plt.figure()
            venn3(sets,mlabels) #DXP=A, GW=B, SP=C
            plt.title(plants[ix]+" "+rank.split("_")[0])
                        
            fig.savefig(plants[ix]+"_"+rank.split("_")[0]+"_venn.png",dpi=1000,bbox_inches="tight")

        
            
            #get shared %
            MP_MG_16S = set.intersection(*sets)
            
            MP_MG =  [i for i in set.intersection( sets[0],sets[1]) if i not in sets[2]]
            MP_16S = [i for i in set.intersection( sets[0],sets[2]) if i not in sets[1]]
            MG_16S = [i for i in set.intersection( sets[1],sets[2]) if i not in sets[0]]        
            
        
        
            contributions=[]
            for file in methodfiles:
                filepath=str(Path(folder,file))
                xlsdf=pd.read_excel(filepath).fillna(0)
                
                for i in [MP_MG_16S,MP_MG,MP_16S,MG_16S]:
                    
                    
                    ratio=sum(xlsdf.loc[xlsdf[rank].isin(i), rank+"_count"])/sum(xlsdf[rank+"_count"])
                    contributions.append(ratio)
                    
            
        
            arr=np.array(contributions).reshape(-1,4)
            df=pd.DataFrame(arr,index=methods,columns=["A_B_C","A_B","A_C","B_C"])
            df["self"]=1-df.sum(axis=1)


    #%%        
    
            colors=np.array([["#d2d2d3","#e7ddcd","#e6cee1","#c6e7e7","#e7b8b9"],
                             ["#d2d2d3","#e7ddcd","#e6cee1","#c6e7e7","#bbdcc2"],
                             ["#d2d2d3","#e7ddcd","#e6cee1","#c6e7e7","#b8c6df"]])
            #sns.set_style("white")
            arr=df.values
            fig, ax = plt.subplots(figsize=[2,4])
            
            xcors=np.linspace(0,0.25*3,3)
            
            for x,i in enumerate(arr):
                bottom=0
                for y,j in enumerate(i):
                    
                    
            
                    ax.bar(xcors[x], height=j, bottom=bottom, width=0.35, color=colors[x,y])
                    #ax.barh(xcors[x], height=j, bottom=bottom, width=0.35, color=colors[x,y])
                    bottom+=j
                    
            #plt.axis('off')
            plt.grid(b=None)
            plt.box(False)
            plt.xticks([])
            ax.set_yticklabels(np.round(ax.get_yticks(),1),rotation=90,va="center")
            ax.yaxis.set_ticks_position("right")
            
            fig.savefig(Path(folder).name+"_"+plants[ix]+"_"+rank.split("_")[0]+"_contribution.png",dpi=1000,bbox_inches="tight")
   
            
        #%%
 
sns.set_style("white")
        

arr=np.array([
["#d2d2d3", "MP U MG U 16S"],
["#e7ddcd", "MP U MG"],
["#e6cee1", "MP U 16S"],
["#c6e7e7", "MG U 16S"],
["#e7b8b9", "MP only"],
["#bbdcc2", "MG only"],
["#b8c6df", "16S only"]])

colors=arr[:,0]
labels=arr[:,1]           

fig, ax = plt.subplots()
for i in colors:
    ax.bar(xcors[x], height=0, bottom=0, width=0, color=i)
plt.legend(labels,loc='center left', bbox_to_anchor=(1, 0.5))
fig.savefig("venn_legend.png",dpi=1000,bbox_inches="tight")
        
        
#matplotlib.pyplot.bar(x, height, width=0.8, bottom=None, *, align='center', data=None, **kwargs)
#stacked bar with matching hex colors
