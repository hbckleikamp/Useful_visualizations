#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 15:47:47 2021

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

#%% Description PCA per plants:
    #-read all plants and techniques, compare
    #-compare abundance profiles separately for different ranks.
    




#%% get Biases

folds=pd.read_excel("/Volumes/Seagate_SSD/third_paper new/Supplemental/3_Difference and fold change/perc_change.xlsx")

#mean perc changes
folds=folds.set_index("taxa")

#16S
cols=[c for c in folds.columns if '16S-MP perc_change' in c]
m=folds[cols].mean(axis=1).fillna(0)
bias_16SMP=m.to_dict()
bias_MP16S=(m*-1).to_dict()

cols=[c for c in folds.columns if '16S-MG perc_change' in c]
m=folds[cols].mean(axis=1).fillna(0)
bias_16SMG=m.to_dict()
bias_MG16S=(m*-1).to_dict()

cols=[c for c in folds.columns if 'MG-MP perc_change' in c]
m=folds[cols].mean(axis=1).fillna(0)
bias_MGMP=m.to_dict()
bias_MPMG=(m*-1).to_dict()


floats=np.array(m.tolist())
colormap = cm.RdYlGn
normalize = mcolors.Normalize(vmin=np.min(floats), vmax=np.max(floats))
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
    


#%% load the 3 MP files and sum them.




class Node:
    def __init__(self,name,rank,count,xcor,ycor):
        self.name=name
        self.rank=rank
        self.count=count
        self.xcor=xcor
        self.ycor=ycor
        #not each node has a parent
        
        
        
folder="/Volumes/Seagate_SSD/third_paper/outputs for figures/branch_quant/"
ranks=["superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]  


plants=["DXP","GW","SP"]
methods=["MPDB","MG","S16"]

biases=[[bias_MPMG,bias_MP16S], #MPDB
        [bias_MGMP,bias_MG16S], #MG
        [bias_16SMP,bias_16SMG] #S16
        ]

#make colorscale


for im,method in enumerate(methods):


    for bias in biases[im]:
        
        
    
    
        
        
        datasets=[]   
        combinedMP=pd.DataFrame()          
        for file in os.listdir(folder):
            if file[0].isalnum() and file.endswith(".xlsx") and method in file:
                
                
        
                plant=[i for i in plants if i in file][0]
                method=[i for i in methods if i in file][0]
                
                vals=pd.read_excel(str(Path(folder,file)))
                vals.pop("Unnamed: 0")
                vals=vals.fillna("")
                combinedMP=pd.concat((combinedMP,vals))
         
        combinedMP=combinedMP.groupby(ranks)["Count"].apply(sum).reset_index()
        
        
        
        def mkempt(x,bias):
            if not x in bias.keys():
                x=""
            return x
            
        combinedMP[ranks]=combinedMP[ranks].applymap(lambda x: mkempt(x,bias))
        
        vals=combinedMP
        emptyinds=np.where(vals.to_numpy()=="")
        for c,i in enumerate(emptyinds[0]):
            vals.iloc[emptyinds[0][c],emptyinds[1][c]]="empty_"+str(c)
            
        
        
        
        #create nodes
        nodes=[]
        
        for rank_ind,rank in enumerate(ranks):
            
            #assign colors based on superkingdom ranks
            color_rank=1
            if rank_ind==color_rank:
                colordict={}
                
                sr_rank_names=vals[ranks[color_rank]].unique().tolist()
                cmap = plt.get_cmap("tab10") #cycle = plt.rcParams['axes.prop_cycle'].by_key()['color'] 
                [colordict.update({i:cmap(c)}) for c,i in enumerate(sr_rank_names)]
            
            counts=vals.groupby(rank)["Count"].sum()/vals["Count"].sum()*100
            names=vals[rank].drop_duplicates()
            ycors=np.linspace(0,10,len(names))
        
            
        
        
            for name_ind,name in enumerate(names):
                if name!="":
                    
                    
                    count=counts.loc[name]
                    if "empty_" in name: count=0
                    
                    
                    mynode=Node(name,rank,count,
                                xcor=rank_ind,
                                ycor=ycors[name_ind])
                    
                    #add parents
                    if rank_ind: 
                        mynode.parent=vals.loc[vals[ranks[rank_ind]]==name,ranks[rank_ind-1]].values[0]
                    
                    #add children
                    if rank_ind<len(ranks)-1: 
                        children=np.unique(vals.loc[vals[ranks[rank_ind]]==name,ranks[rank_ind+1]].values).tolist()
                
                        mynode.children=[i for i in children if i !=""]
                    
        
                        
                    
                    
                    nodes.append(mynode)
                      
                           
              
               
        # center parents to children
        for rank in ranks[::-1]:
            
            for node in nodes:
                if node.rank==rank:
                    if hasattr(node, 'children'):
                        if len(node.children):
           
                            ycors=[]
                            for child in node.children:
        
                                ycors.append([node.ycor for node in nodes if node.name==child][0])
                            
                            node.ycor=(max(ycors)+min(ycors))/2
                            node.childycors=ycors
                            
                            
        
                            
                
        
        
        
        # add parent ycor               
        
        for node in nodes:
            if hasattr(node, 'parent'):
                parent_name=node.parent
                parentycor=[node.ycor for node in nodes if node.name==parent_name]
                node.parentycor=parentycor
              
        
        # add colors
        for node in nodes:
            if bias.get(node.name):
                node.color=s_map.to_rgba(bias.get(node.name))
            else:
                node.color=s_map.to_rgba(0)
          
        
        #scale figsize (ylen) with len of vals?
        figsize=(60,40)  #xlen, ylen
        plt.figure(figsize=figsize) #scale thickness with figsize
        for node in nodes:
         
            
            if hasattr(node, 'parent'):
                
                point1=[node.xcor,node.ycor]#child
                point2=[node.xcor-1,node.parentycor] #parent
        
                #plot hanging line
                a = (point2[1] - point1[1])/(np.cosh(point2[0]) - np.cosh(point1[0]))
                b = point1[1] - a*np.cosh(point1[0])
                x = np.linspace(point1[0], point2[0], 100)
                y = a*np.cosh(x) + b
                
                plt.plot(x, y, c=node.color,linewidth=node.count*0.45*(figsize[1]/4), alpha=0.8)
        
                
        
        for node in nodes:
            plt.plot(node.xcor,node.ycor,".",color=node.color,markersize=node.count*figsize[1]/4, alpha=0.8)
        
            if "empty_" not in node.name and node.count>0.5:      
                plt.text(node.xcor, node.ycor, node.name, fontsize=2*(figsize[1]/6))
                plt.text(node.xcor-0.1, node.ycor, round(node.count,2), fontsize=2*(figsize[1]/6))
                
            
         
        plt.xlim([-0.52,6.2])
        
        plt.axis('off')
        plt.grid(b=None)
        

        plt.title(plant+method)         
        plt.show()


           
