#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 13:49:00 2022

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


gdfs=pd.read_csv("/Volumes/Seagate_SSD/paper 3/third_paper_newest/Figures_main/4 Genus Metabolic /found_proteins_integrated_annotation.tsv",sep="\t")    
#%% difference



#lets properly sort these categories


# J	FCCCFC	Translation, ribosomal structure and biogenesis
# A	FCDCFC	RNA processing and modification
# K	FCDCEC	Transcription
# L	FCDCDC	Replication, recombination and repair
# B	FCDCCC	Chromatin structure and dynamics
# D	FCFCDC	Cell cycle control, cell division, chromosome partitioning
# Y	FCFCCC	Nuclear structure
# V	FCFCBC	Defense mechanisms
# T	FCFCAC	Signal transduction mechanisms
# M	ECFCAC	Cell wall/membrane/envelope biogenesis
# N	DCFCAC	Cell motility
# Z	CCFCAC	Cytoskeleton
# W	BCFCAC	Extracellular structures
# U	ACFCAC	Intracellular trafficking, secretion, and vesicular transport
# O	9CFCAC	Posttranslational modification, protein turnover, chaperones
# X	9CFC9C	Mobilome: prophages, transposons
# C	BCFCFC	Energy production and conversion
# G	CCFCFC	Carbohydrate transport and metabolism
# E	DCFCFC	Amino acid transport and metabolism
# F	DCECFC	Nucleotide transport and metabolism
# H	DCDCFC	Coenzyme transport and metabolism
# I	DCCCFC	Lipid transport and metabolism
# P	CCCCFC	Inorganic ion transport and metabolism
# Q	BCCCFC	Secondary metabolites biosynthesis, transport and catabolism
# R	E0E0E0	General function prediction only
# S	CCCCCC	Function unknown

# cats=[      
# # 'AOB','NOB','DNR','PAO','GAO', #nutrient removal

# "M","N","W","T","U","porin",

# # "Gt","Et","Ft","Ht","It","Pt","Qt", #membrane/extracllular
# # "C","Gm","Em","Fm","Hm","Im","Pm","Qm", #metabolism
# "C","G","E","F","H","I","P","Q",
# "J","A","K","L","O","B","D","Y","Z", #housekeeping
# "V","X" #other
# ]

cats=[      
# 'AOB','NOB','DNR','PAO','GAO', #nutrient removal

"Nm","C","G","E","F","H","I","P","Q",

"M","N","W","T","U","porin",

# "Gt","Et","Ft","Ht","It","Pt","Qt", #membrane/extracllular


"J","A","K","L","O","B","D","Y","Z", #housekeeping
"V","X","R","S" #other
]

#catgory descriptions
metabolism=["Nitrogen", "Energy" ,"Carbohydrate","Amino acid","Nucleotide","Coenzyme","Lipid","Inorganic ion","Secondary metabolite"]
membrane=["Membrane","Motility","Extracellular structures","Signal transduction","Secretion","Porin"]
cell_cycle=["Translation","Transcription","Replication","PTMs & chaperones","Chromatin","Cell cycle","Cytoskeleton"]
other=["Defense","Mobilome","Predictive","Unknown"]
cat_desc=metabolism+membrane+cell_cycle+other


#%%
#fraction based on mean
t=gdfs[["Accession",'PSMs','PSMs_normalized','AvgDepth','Depth_normalized', 'PSM/Depth',"cog_cat"]].drop_duplicates()
t["AvgDepth"]=t["AvgDepth"].replace("",0).astype(float)


r=t.groupby("cog_cat")[['PSMs','AvgDepth']].sum()
r=r[r.index!=""]

spsm=r["PSMs"].sort_values()
mpsm=spsm.median()
rpsm=spsm/mpsm


sd=r['AvgDepth'].sort_values()
md=sd.median()
rd=sd/md        #median ratio
medrad=rpsm/rd
medrad=medrad.sort_values(ascending=False)



#%%
#common fraction
t=gdfs[["Accession",'PSMs','PSMs_normalized','Depth_normalized', 'PSM/Depth',"cog_cat"]].drop_duplicates()
r=t.groupby("cog_cat")[['PSMs_normalized','Depth_normalized']].sum()
r=r[r.index!=""]
r=r/r.sum()*100
cats=[c for c in cats if c in r.index]
#diff=r["PSMs_normalized"]-r["Depth_normalized"]
#diff=(r["PSMs_normalized"]-r["Depth_normalized"])/((r["PSMs_normalized"]+r["Depth_normalized"])/2)
diff=(rpsm-rd)/((rpsm+rd)/2) #percentage difference from fold change of the median
diff=diff.loc[cats]


#rename cats
diff=pd.DataFrame(diff,columns=["diff"])
diff["desc_cat"]=diff.index+[": "]*len(cat_desc)+cat_desc
diff["desc"]=cat_desc

#group cats
ds=[]
for i in [metabolism,membrane,cell_cycle,other]:
    
    ds.append(diff[diff["desc"].isin(i)].sort_values(by="diff",ascending=False))
    
ds=pd.concat(ds)
    
#sort by group
#ds["desc"]=ds.index+[": "]*len(cat_desc)+cat_desc
ds=ds.set_index("desc_cat")
ds=pd.Series(ds["diff"])

#add back column descriptors

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm


ind=ds.index
dim=(len(ind),10)
names=ds.index[::-1]
floats=ds.values
colormap = cm.RdYlGn
normalize = mcolors.Normalize(vmin=np.min(floats), vmax=np.max(floats))
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)

fig = plt.figure()
defsize = fig.get_size_inches()*fig.dpi
sns.set_style("white")
fig, ax = plt.subplots(figsize=(dim[1]*0.5,dim[0]/3))
width=0.45
bottom=0
for x, i in enumerate(ds[::-1]):
    ax.barh(x, i, width, bottom, color=s_map.to_rgba(i))#, hatch=patterns[0],zorder=5)#, yerr=std_allMP16S[x] )

ax.set_yticks(np.arange(len(names)))
ax.set_yticklabels(names, Fontsize=9, rotation=-40,rotation_mode="anchor",verticalalignment="top")
ax.xaxis.tick_top()
ax.yaxis.tick_left()
ax.tick_params(axis='y',length=2)

ax.set_ylim([-0.5,len(ds)-0.2])
fig.savefig("common_cog_cats_difference.png",dpi=400,bbox_inches="tight")


#%%