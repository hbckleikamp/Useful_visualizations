#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:18:37 2021

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
 
sns.set_style("white")
       

folder="/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/main_ouput/"
ranks=["superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]  
plants=["DXP","GW","SP"]
methods=["MP","MG","S16"]


class datapoint:
        def __init__(self,name,rank,plant,method,richness,evenness,diversity):
                        
            self.name=name
            self.rank=rank
            self.plant=plant
            self.method=method
            self.richness=richness
            self.evenness=evenness
            self.diversity=diversity





import skbio

datapoints=[]
for method in methods:
    

    
    
    for file in os.listdir(folder):
        if file[0].isalnum() and file.endswith(".xlsx") and method in file:
            
            data=pd.read_excel(str(Path(folder,file)))
            plant=[i for i in plants if i in file][0]
            method=[i for i in methods if i in file][0]
            
    
            for rank in ranks:
                
                name=plant+"_"+method+"_"+rank
                dat=data[rank+"_count"].dropna().tolist()
                
                richness=sum(pd.notnull(data[rank]))
                evenness=skbio.diversity.alpha.simpson_e(dat) #simpson eveness 0-1
                diversity=skbio.diversity.alpha_diversity('shannon',dat)[0]
                
                datapoints.append(datapoint(name,rank,plant,method,richness,evenness,diversity))


#%%

names=[]

avrich=[]
stdrich=[]

aveven=[]
stdeven=[]

avdivs=[]
stddivs=[]

ranks=["phylum_name","class_name","order_name","family_name","genus_name","species_name"]  
for rank in ranks:
    

        
    
    
    for method in methods:
        rich=[]
        even=[]
        div=[]
        
        
        for plant in plants:
            
            
            
            
    
            for datapoint in datapoints:
                if datapoint.plant==plant and datapoint.method==method and datapoint.rank==rank: #genus only?
                
                    
                    names.append(datapoint.name)
                    rich.append(datapoint.richness)
                    even.append(datapoint.evenness)
                    div.append(datapoint.diversity)
                
        
            
            
        avrich.append(np.mean(rich))
        stdrich.append(np.std(rich))
        
        aveven.append(np.mean(even))
        stdeven.append(np.std(even))
        
        avdivs.append(np.mean(div))
        stddivs.append(np.std(div))
            
            
            
    cmap=plt.cm.get_cmap('Pastel1')#'tab10')
    colors=[]
    for method in methods:
        # if "S16" in method: colors.append("#b8c6df")#cmap(2))
        # if "MP" in method: colors.append("#e7b8b9")#cmap(1))
        # if "MG" in method: colors.append("#bbdcc2")#cmap(4))
        if "S16" in method: colors.append(cmap(1))
        if "MP" in method: colors.append(cmap(0))
        if "MG" in method: colors.append(cmap(2))

    

        
def grouped_bar(labels,title,ylabel,groups_of,data,err,color=colors):
    rs=np.array(data).reshape(-1,groups_of)
    es=np.array(err).reshape(-1,groups_of)
    
    x = np.arange(len(labels))
    width = 1/groups_of*0.6
    
    fig, ax = plt.subplots()
    
    xrange=np.linspace(x - width, x + width, groups_of).T
    
    
    for i in range(len(labels)):
        ax.bar(xrange[i], rs[i], width, color=colors,
               yerr=es[i])
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    fig.suptitle(title)



labels=[rank.split("_")[0] for rank in ranks]

  
    
#%%

rs=np.array(avrich).reshape(-1,3).T
es=np.array(stdrich).reshape(-1,3).T

fig, ax = plt.subplots()
for ix,i in enumerate(rs):
    plt.errorbar(np.arange(len(rs[ix])), rs[ix], yerr=es[ix], color=colors[ix])
    
    data = {
        'x': np.arange(len(rs[ix])),
        'y1': [y - e for y, e in zip(rs[ix], es[ix])],
        'y2': [y + e for y, e in zip(rs[ix], es[ix])]}
    plt.fill_between(**data, 
                     alpha=.40, 
                     color=colors[ix])
    
ylabel="number of taxa"
ax.set_ylabel(ylabel, fontsize=12)
ax.set_xticks(np.arange(len(rs[ix])))
ax.set_xticklabels(labels, rotation=40, fontsize=12)

plt.title('Richness')
fig.savefig("Richness.png",dpi=1000,bbox_inches="tight")

rs=np.array(aveven).reshape(-1,3).T
es=np.array(stdeven).reshape(-1,3).T

fig, ax = plt.subplots()
for ix,i in enumerate(rs):
    plt.errorbar(np.arange(len(rs[ix])), rs[ix], yerr=es[ix], color=colors[ix])
       
    data = {
        'x': np.arange(len(rs[ix])),
        'y1': [y - e for y, e in zip(rs[ix], es[ix])],
        'y2': [y + e for y, e in zip(rs[ix], es[ix])]}
    plt.fill_between(**data, 
                     alpha=.40, 
                     color=colors[ix])
    
ylabel="Simpson's evenness"
ax.set_ylabel(ylabel, fontsize=12)
ax.set_xticks(np.arange(len(rs[ix])))
ax.set_xticklabels(labels, rotation=40, fontsize=12)
    
plt.title('Evenness')
fig.savefig("Simpson's evenness.png",dpi=1000,bbox_inches="tight")

rs=np.array(avdivs).reshape(-1,3).T
es=np.array(stddivs).reshape(-1,3).T

fig, ax = plt.subplots()
for ix,i in enumerate(rs):
    plt.errorbar(np.arange(len(rs[ix])), rs[ix], yerr=es[ix], color=colors[ix])
    
    data = {
        'x': np.arange(len(rs[ix])),
        'y1': [y - e for y, e in zip(rs[ix], es[ix])],
        'y2': [y + e for y, e in zip(rs[ix], es[ix])]}
    plt.fill_between(**data,
                     alpha=0.4#.20
                     , color=colors[ix])
    
    
ylabel="Shannon diversity"
ax.set_ylabel(ylabel, fontsize=12)
ax.set_xticks(np.arange(len(rs[ix])))
ax.set_xticklabels(labels, rotation=40, fontsize=12)

plt.title('Diversity')
fig.savefig("Shannon diversity.png",dpi=1000,bbox_inches="tight")

#%% legend plot

sns.set_style("white")

# if "S16" in method: colors.append(cmap(1))
# if "MP" in method: colors.append(cmap(0))
# if "MG" in method: colors.append(cmap(2))
        

arr=np.array([
[cmap(0), "MP"],
[cmap(2), "MG"],
[cmap(1), "16S"]
])

colors=arr[:,0]
labels=arr[:,1]           

fig, ax = plt.subplots()
for i in colors:
    ax.bar([0,0], height=0, bottom=0, width=0, color=i)
plt.legend(labels,loc='center left', bbox_to_anchor=(1, 0.5))
fig.savefig("alpha_legend.png",dpi=1000,bbox_inches="tight")
        