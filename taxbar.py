#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 22:41:28 2021

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
# sns.set_style("white")

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
print(os.getcwd())

#%% grouped bar per rank


import pandas as pd
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt

import itertools
def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])





#%% PCA (different sizes)
ranks=["superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]  
plants=["DXP","GW","SP"]
methods=["MPDB","MG","S16"]


folder="/Volumes/Seagate_SSD/third_paper new/Annotation/Output files/rank_quant/"

for rank in ranks[1:]:#ranks[1:-1]:
        
    once=True 
    for method in methods:
        for plant in plants:
            for file in os.listdir(folder):
                if file[0].isalnum() and file.endswith(".xlsx") and method in file and plant in file:
    
                    
                    pdf=pd.read_excel(str(Path(folder,file)))[[rank,rank+"_count"]].dropna()
                    pdf.columns=["name","count "+plant+" "+method]
                    
                    
                    
                    if once:
                        df=pdf
                        once=False
                    else:
                        df=df.merge(pdf,on="name",how="outer")
                    
        
            
    df["name"]=df["name"].apply(lambda x: x[3:])
    df=df.set_index("name").fillna(0)
    
    top=df.mean(axis=1).sort_values(ascending=False).index[0:8]
    tdf=df[df.index.isin(top)].append(pd.DataFrame(df[~df.index.isin(top)].sum(),columns=["other"]).T)
    
    #first plots
    # figure, axes = plt.subplots()
    # ax=tdf.T.plot.bar(stacked=True,legend=False)
    
    
    
    # figure, axes = plt.subplots()
    # ax=(tdf.T*0.01).plot.bar(stacked=True,legend=False)
    # ax.set_xlim([10,100])
    # ax.legend(tdf.index,facecolor='white',loc="right")
    
    
    #https://stackoverflow.com/questions/22787209/how-to-have-clusters-of-stacked-bars-with-python-pandas
    # create fake dataframes
    df1 = pd.DataFrame(tdf[tdf.columns[0:len(tdf.columns):3]]).T
    df2 = pd.DataFrame(tdf[tdf.columns[1:len(tdf.columns):3]]).T
    df3 = pd.DataFrame(tdf[tdf.columns[2:len(tdf.columns):3]]).T
    
    df1.index=["MP","MG","16S"]
    df2.index=["MP","MG","16S"]
    df3.index=["MP","MG","16S"]

    
    
    dfall=(df1,df2,df3)
    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    
    figure, axes = plt.subplots()
    axe = plt.subplot(111)
    
    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      ax=axe,
                      legend=False) # make bar plots
    
    
    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))   
                rect.set_width(1 / float(n_df + 1)*0.9)
    
    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 0)
    
    
    l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0.2])
    axe.set_title(rank)
    handles, labels = axe.get_legend_handles_labels()
    plt.legend(flip(h[:n_col], 2), flip(l[:n_col], 2), loc=[-0.1, -0.4], ncol=3)

        
