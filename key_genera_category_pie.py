#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 14:31:01 2021

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

from io import StringIO

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
print(os.getcwd())


#%% key genera

pfiles=[
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/Proteomics/db_comp_rank_quantf0.1/DXP_MP_quant_total_CBFB_5_0.5_contig_lca_proteins_blca_GTDB16S_all_FL_hom_tax_merged_proteins.tsv.xlsx",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/Proteomics/db_comp_rank_quantf0.1/GW_MP_quant_total_CBFB_5_0.5_contig_lca_proteins_blca_GTDB16S_all_FL_hom_tax_merged_proteins.tsv.xlsx",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/Proteomics/db_comp_rank_quantf0.1/SP_MP_quant_total_CBFB_5_0.5_contig_lca_proteins_blca_GTDB16S_all_FL_hom_tax_merged_proteins.tsv.xlsx"
]
        


ranks=["superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]   
samples=["DXP","GW","SP"]

ugs=[]
for ix,p in enumerate(pfiles):
    

    df=pd.read_excel(p)
    ug=pd.DataFrame(df[df["genus_name_count"]>3]["genus_name"])
    ug["plant"]=samples[ix]
    

    ugs.append(ug)
    

ugs=pd.concat(ugs)
ugs.index=ugs["genus_name"]+" "+ugs["plant"]

#%%Sort by bias
change="/Volumes/Seagate_SSD/paper 3/third_paper_newest/Figures_main/4 Genus Metabolic /perc_change.xlsx"
cdf=pd.read_excel(change)

cdf=cdf[cdf["taxa"].str.startswith("g__")]
cdf=cdf.set_index("taxa")


#shared absolute difference
once3=True
for sample in samples:
    sh=cdf[[sample +" MG",sample +" MP"]]
    sh=sh[(sh>0).all(axis=1)]
    sh=sh/sh.sum()*100

    if once3:
        shs=sh
        once3=False
    else:
        shs=pd.concat([shs,sh],axis=1)

shs=shs.fillna(0)

diffs=[]
for ix in ugs.index:
    
    sample=ix.split(" ")[-1]
    gen=ix.replace(" "+sample,"")
    diffs.append(shs.loc[gen,sample+ " MG"]-shs.loc[gen,sample+ " MP"])
    
ugs["diffs"]=np.array(diffs)



#%% cog cats

key_g=pd.read_csv("/Volumes/Seagate_SSD/paper 3/third_paper_newest/Figures_main/4 Genus Metabolic /key genera.csv")
gdfs=pd.read_csv("/Volumes/Seagate_SSD/paper 3/third_paper_newest/Figures_main/4 Genus Metabolic /found_proteins_integrated_annotation.tsv",sep="\t")

cats=[      "Nm","C","G","E","F","H","I","P","Q",
"M","N","W","T","U","porin",
"J","K","L","O","B","D","Z", #housekeeping
"V","X","R","S" #other
]

#catgory descriptions
metabolism=["Nitrogen", "Energy" ,"Carbohydrate","Amino acid","Nucleotide","Coenzyme","Lipid","Inorganic ion","Secondary metabolite"]
membrane=["Membrane","Motility","Extracellular structures","Signal transduction","Secretion","Porin"]
cell_cycle=["Translation","Transcription","Replication","PTMs & chaperones","Chromatin","Cell cycle","Cytoskeleton"]
other=["Defense","Mobilome","Predictive","Unknown"]
cat_desc=metabolism+membrane+cell_cycle+other

df=pd.DataFrame(cats)
df.columns=["cats"]
df["desc_cat"]=df["cats"]+[": "]*len(cat_desc)+cat_desc
df=df.set_index("cats")

#%%
#fraction based on mean
t=gdfs[["Accession",'PSMs','PSMs_normalized','AvgDepth','Depth_normalized', 'PSM/Depth',"cog_cat"]].drop_duplicates()
r=t.groupby("cog_cat")[['PSMs','AvgDepth']].sum()
r=r[r.index!=""]

spsm=r["PSMs"].sort_values()
mpsm=spsm.median()
rpsm=spsm/mpsm

sd=r['AvgDepth'].sort_values()
md=sd.median()
rd=sd/md

medrad=rpsm/rd
medrad=medrad.sort_values(ascending=False)



#%%
#common fraction
t=gdfs[["Accession",'PSMs','PSMs_normalized','Depth_normalized', 'PSM/Depth',"cog_cat"]].drop_duplicates()
d=t.groupby("cog_cat")[['PSMs_normalized','Depth_normalized']].sum()
d=d[d.index!=""]
d=d/d.sum()*100


#top cats

ind=d.sum(axis=1).sort_values(ascending=False).index
top10=ind[0:10].tolist()
oind=ind[10:].tolist()

other=d.loc[ind[10:]]


d=d.loc[top10]

d["desc_cat"]=df.loc[d.index,"desc_cat"]


a=other.sum().values.tolist()+["other"]
d.loc["other",:]=a
cmap = plt.get_cmap("tab10")
d["colors"]=[cmap(i) for i in range(10)]+["thistle"]


#%%
# key genera

pfiles=[
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/Proteomics/db_comp_rank_quantf0.1/DXP_MP_quant_total_CBFB_5_0.5_contig_lca_proteins_blca_GTDB16S_all_FL_hom_tax_merged_proteins.tsv.xlsx",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/Proteomics/db_comp_rank_quantf0.1/GW_MP_quant_total_CBFB_5_0.5_contig_lca_proteins_blca_GTDB16S_all_FL_hom_tax_merged_proteins.tsv.xlsx",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/Proteomics/db_comp_rank_quantf0.1/SP_MP_quant_total_CBFB_5_0.5_contig_lca_proteins_blca_GTDB16S_all_FL_hom_tax_merged_proteins.tsv.xlsx"
]
        


ranks=["superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]   
samples=["DXP","GW","SP"]

ug=key_g["taxa"]
# ug=["g__Nitrosomonas","g__Azonexus","g__Accumulibacter","g__Competibacter_A","g__Nitrospira_A","g__propionivibrio"]

groups=ugs.groupby("plant")


PSMs=[]
Depths=[]

for g in groups:
    
    g=g[1]
    
    genera=g["genus_name"]
    plant=g["plant"].iloc[0]
    
    for genus in genera:
       
        
        
        t=gdfs.loc[(gdfs["Plant"]==plant) & (gdfs["genus_name"]==genus),
                    ["Accession",'PSMs','PSMs_normalized','Depth_normalized', 'PSM/Depth',"cog_cat"]].drop_duplicates()
        
        if len(t):
            r=t.groupby("cog_cat")[['PSMs_normalized','Depth_normalized']].sum()
            r=r[r.index!=""]
            r=r/r.sum()*100
    
    
            o=r[r.index.isin(oind)]
            r=r[r.index.isin(top10)]
            r.loc["other",:]=o.sum()
            
            PSMs.append(r["PSMs_normalized"])
            Depths.append(r["Depth_normalized"])
            
            fig = plt.figure()
            r["PSMs_normalized"].plot.pie(colors=d.loc[r.index,"colors"].tolist())
            plt.title(plant+"_"+genus+" PSMs")
            fig.savefig(plant+"_"+genus+"_cog_cats_PSMs_Pie.png",dpi=400,bbox_inches="tight")
        
            fig = plt.figure()
            r["Depth_normalized"].plot.pie(colors=d.loc[r.index,"colors"].tolist())
            plt.title(plant+"_"+genus+" Depth")
            fig.savefig(plant+"_"+genus+"_cog_cats_Depth_Pie.png",dpi=400,bbox_inches="tight") 
    
            #slim output for custom mergeing
            
            fig = plt.figure()
            r["PSMs_normalized"].plot.pie(colors=d.loc[r.index,"colors"].tolist(),title='',ylabel='',fontsize=0)
            fig.savefig(plant+"_"+genus+"_cog_cats_PSMs_Pie_slim.png",dpi=400,bbox_inches="tight")

        
            fig = plt.figure()
            r["Depth_normalized"].plot.pie(colors=d.loc[r.index,"colors"].tolist(),title='',ylabel='',fontsize=0)
            fig.savefig(plant+"_"+genus+"_cog_cats_Depth_Pie_slim.png",dpi=400,bbox_inches="tight")


 #%%
sns.set_style("white")
labels=ugs["genus_name"]+" "+ugs["plant"]
PSMbar=pd.concat(PSMs,axis=1).fillna(0).T
Dbar=pd.concat(Depths,axis=1).fillna(0).T
PSMbar.index=labels
Dbar.index=labels

PSMbar=PSMbar.loc[labels.sort_values()]
Dbar=Dbar.loc[labels.sort_values()]

f = plt.figure()

plt.title('PSMs expression', color='black')
PSMbar.plot(kind='bar', stacked=True, ax=f.gca(),color=d.loc[r.index,"colors"])
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel("COG category PSM distribution")
f.savefig("PSMs_expression.png",dpi=1000,bbox_inches="tight")

f = plt.figure()

plt.title('Depth expression', color='black')
Dbar.plot(kind='bar', stacked=True, ax=f.gca(),color=d.loc[r.index,"colors"])
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel("COG category read distribution")
f.savefig("Depth_expression.png",dpi=1000,bbox_inches="tight")

#%%

f=plt.figure()
diffs=pd.DataFrame(ugs.loc[labels.sort_values(),"diffs"]).T
diffs.columns=np.round(diffs.iloc[0,:].values,1).tolist()
g=sns.heatmap(diffs,cmap='RdYlGn',center=0,annot=True)
g.set_xticklabels(g.get_xticklabels(), 
                  rotation=0,fontsize=len(diffs.iloc[0])/4)
f.savefig("expression_bias.png",dpi=1000,bbox_inches="tight")

#%% plot legend


sns.set_style("white")
colors=d["colors"].tolist()
fig = plt.figure(figsize=None)
for c in colors:
    plt.bar(0,0,color=c)


plt.legend(d["desc_cat"].tolist())
fig.savefig('status_legend.png',dpi=400,bbox_inches="tight")