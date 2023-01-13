#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 20:39:38 2021

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


ranks=["superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]
plants=["DXP","GW","SP"]

once=True


#%% get genus lvl hits

MPDFB_files=[
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/Proteomics/db_comp_rank_quantf0.1/DXP_MP_quant_total_CBFB_5_0.5_contig_lca_proteins_blca_GTDB16S_all_FL_hom_tax_merged_proteins.tsv.xlsx",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/Proteomics/db_comp_rank_quantf0.1/GW_MP_quant_total_CBFB_5_0.5_contig_lca_proteins_blca_GTDB16S_all_FL_hom_tax_merged_proteins.tsv.xlsx",
"/Volumes/Seagate_SSD/paper 3/third_paper_newest/Annotation_Quantification/2_Quantification/Proteomics/db_comp_rank_quantf0.1/SP_MP_quant_total_CBFB_5_0.5_contig_lca_proteins_blca_GTDB16S_all_FL_hom_tax_merged_proteins.tsv.xlsx"


    ]



c="/Volumes/Seagate_SSD/paper 3/third_paper_newest/Figures_main/4 Genus Metabolic /diff.xlsx"


#%%
f="/Volumes/Seagate_SSD/paper 3/third_paper_newest/Figures_main/4 Genus Metabolic /found_proteins_integrated_annotation.tsv"
gdfs=pd.read_csv(f,sep="\t")
cdf=pd.read_excel(c).fillna(0)
cdf=cdf[cdf["taxa"].str.startswith("g__")]
cdf=cdf.set_index("taxa")
#%%


#%%
def lca(lin,ranks=ranks):
            #lca
        lca=[]
        #print(lin)
        if len(lin)>1:
            #lin=pd.DataFrame(lin,columns=ranks)
            
            ix=0
            for r in ranks:
                if lin[r].nunique()!=1:
                    break
                else:
                    ix+=1
            lca=lin.iloc[0,0:ix].tolist()
        elif len(lin)==1:
            lca=lin.iloc[0].tolist() #since this produces a list
        
        lca=lca+[" "]*(len(ranks)-len(lca))
   
            
        return lca


#perform "lca"

lcadf=gdfs.groupby("Accession")[ranks].apply(lambda x: lca(x))
lcadf=pd.DataFrame(lcadf.tolist(),index=lcadf.index).reset_index()
lcadf.columns=["Accession"]+ranks
fundf=gdfs[["Accession","PSMs","role_genes","Plant"]].drop_duplicates().merge(lcadf,on='Accession',how='left')                                             
#%%
#only select genera that show up in two techniques.
genera=[]
for plant in plants:
    pc=[c for c in cdf.columns if plant in c][0:3]
    tx=cdf.loc[(cdf[pc]>0).sum(axis=1)>0].index
    genera.extend(tx.tolist())
genera=list(set(genera))




#%% establish fold change order

#shared absolute difference
once3=True
for plant in plants:
    sh=cdf[[plant +" MG",plant +" MP"]]
    sh=sh[(sh>0).all(axis=1)]
    sh=sh/sh.sum()*100
    shdiff=sh[plant +" MG"]-sh[plant +" MP"]

    if once3:
        sdifdf=shdiff
        shs=sh
        once3=False
    else:
        sdifdf=pd.concat([sdifdf,shdiff],axis=1)
        shs=pd.concat([shs,sh],axis=1)
    

MPMGdiff=pd.DataFrame(sdifdf.mean(axis=1),columns=["MPMG"])


#shared absolute difference
once3=True
for plant in plants:
    sh=cdf[[plant +" S16",plant +" MP"]]
    sh=sh[(sh>0).all(axis=1)]
    sh=sh/sh.sum()*100
    shdiff=sh[plant +" S16"]-sh[plant +" MP"]

    if once3:
        sdifdf=shdiff
        shs=sh
        once3=False
    else:
        sdifdf=pd.concat([sdifdf,shdiff],axis=1)
        shs=pd.concat([shs,sh],axis=1)
    

MP16Sdiff=pd.DataFrame(sdifdf.mean(axis=1),columns=["MP16S"])

both=pd.concat([MPMGdiff,MP16Sdiff],axis=1)
#sort_ind=both.mean(axis=1).sort_values().index
sort_ind=MPMGdiff.sort_values(by='MPMG').index #sort only on MGMP diff

both=both.loc[sort_ind]

allMPMG=both.loc[sort_ind]["MPMG"].fillna(0)
allMP16S=both.loc[sort_ind]["MP16S"].fillna(0)



#%% get genus abundance

once=True
for sample in ["MP","MG","S16"]:
    cols=[m+ " "+ sample for m in ["DXP","GW","SP"]]
    dat=cdf.loc[:,cols]
    dat=dat.fillna(0)
    
    if once:
        
        allcomb=dat
        once=False
    else:
        allcomb=pd.concat([allcomb, dat], axis=1)
        
allcomb=allcomb.fillna(0)
allcomb=allcomb[allcomb.index.isin(genera)]
allcomb=allcomb/allcomb.sum()*100

#allowed genera:
    #-expressed nutrient removal genes/or has more than 3% abundance in one measurement
    #-detected  oncein at least two technques
allowed_genera=list(set(fundf.loc[fundf['role_genes'].notnull(),"genus_name"].dropna().tolist()
                        +allcomb[allcomb.max(axis=1)>3].index.tolist()))

sgenus=allcomb.index[allcomb.index.isin(allowed_genera)]
sgenus=both[both.index.isin(sgenus)].index
allcomb=allcomb.loc[sgenus] #sorted

#allcomb=allcomb.loc[both.index]

# reshape to 3s

def chunks(lst,n):
    for i in range(0,len(lst),n):
        yield lst[i:i+n]
        
cols=chunks(allcomb.columns,3)


methods=["MP","MG","S16"]
plants=["DXP","GW","SP"]
once=True
for ix,col in enumerate(cols):
    rdf=pd.DataFrame(allcomb[col].values.flatten(),columns=[methods[ix]])
    rind=[]
    [[rind.extend(i+" "+p for p in plants)] for i in allcomb.index]
    rdf.index=rind
    if once:
        
        rallcomb=rdf
        once=False
    else:
        rallcomb=pd.concat([rallcomb, rdf], axis=1)


sns.set_style("white")
fig, ax = plt.subplots(figsize=(4, 12))
yticklabels=3
xticklabels=3
g=sns.heatmap(rallcomb,
               #cmap="YlGnBu",
               center=20,
               #cmap="vlag",
               mask=(rallcomb==0),
               yticklabels=yticklabels,
               robust=True,
               vmin=0, vmax=40,
               cbar='right'
               )




g.set_yticklabels([str(i).split("'")[1].split(" ")[0].split("__")[1] for i in g.get_yticklabels()],
                  fontsize=(fig.get_size_inches()*fig.dpi)[0]/19,
                  va='top')
g.set_xticklabels(g.get_xticklabels(), 
                  #rotation=30, 
                  fontsize = (fig.get_size_inches()*fig.dpi)[1]/50)


#add grid
figs=(fig.get_size_inches()*fig.dpi).prod() #pixel area
h = g.get_yticks()
h=np.array([0]+h.tolist()+[max(h)+yticklabels]) #add zero line
w = g.get_xticks()
w=np.array([0]+w.tolist()+[max(w)+xticklabels]) #add zero line

#adjust cutoffs manually to fit grid
g.hlines(h-0.6, w.max(),w.min(), linestyle='-', linewidth=1, color="black")
g.vlines(w-0.51,h.max(),h.min(), linestyle='-', linewidth=1, color="black")

plt.tick_params(axis='both', which='major', labelbottom = False, bottom=False, top = False, labeltop=True)
fig.savefig("key genera.png",dpi=400,bbox_inches="tight")
allcomb.to_csv("key genera.csv")
#plt.title(lca+" Abundance Difference")
#plt.xlabel("absolute abundance difference" )
#plt.ylabel("genera" )
#fig.savefig(lca+"_genera.png",dpi=400,bbox_inches="tight")

#%% plot genes


exp=fundf.groupby(["genus_name","Plant","role_genes"])["PSMs"].sum().reset_index()


#normalize per plant
for plant in plants:
    exp.loc[exp["Plant"]==plant,"PSMs"]=exp.loc[exp["Plant"]==plant,"PSMs"]/exp.loc[exp["Plant"]==plant,"PSMs"].sum()


gene_names=["ppk","ppa","blg","glg","Hao","Amo","Nxr","NirS_NirK","Nor","Nos","Nif","Nar_Nap","Nir","Nrf","hzs","hdh","Cyc"][::-1]
             #PAO          GAO          AOB     NOB    "Denitrification"      "Fixation" "Ammonification" "Annamox" 
  
#gene_names=["Hao","Amo","Nxr","Nar_Nap","NirS_NirK","Nor","Nos","Nrf_Nir_Cyc","ppk","ppa","Pho","blg","glg","hzs","hdh"]
#genes=[Hao,Amo,Nxr,Nar_Nap,NirS_NirK,Nos,Nor_Nrf_Nir_Cyc,ppk,ppa,Pho,blg,glg,hzs,hdh]

# gene_names=["ppk","ppa","Pho","blg","glg","Hao","Amo","Nxr","Nar_Nap","NirS_NirK","Nor","Nos","Nor_Nrf_Nir_Cyc","hzs","hdh"][::-1]


gdf=pd.DataFrame(gene_names)
gdf=gdf[gdf.isin(exp["role_genes"].tolist())].dropna()
gdf.columns=["role_genes"]



for genus in sgenus:
    
    g=exp[exp["genus_name"]==genus]
    
    for plant in plants:

        p=g[g["Plant"]==plant][["role_genes","PSMs"]]
        p=p.set_index("role_genes")
        p.columns=[genus +" "+ plant]
        gdf=gdf.merge(p,on="role_genes",how="left")

gdf=gdf.fillna(0)        
       
gdf=gdf.set_index("role_genes")


sns.set_style("white")
fig, ax = plt.subplots(figsize=(16, 12))
yticklabels=1
xticklabels=3
g=sns.heatmap(gdf,
               #cmap="BuPu",
               mask=(gdf==0),
               yticklabels=yticklabels,
               xticklabels=xticklabels,
               center=0.1,
               robust=True,
               cbar=None,
               )


g.set_xticklabels([str(i).split("'")[1].split(" ")[0].split("__")[1] for i in g.get_xticklabels()],
                  fontsize=(fig.get_size_inches()*fig.dpi)[0]/50,
                  va='top')
g.set_yticklabels(g.get_yticklabels(), rotation=60, fontsize = (fig.get_size_inches()*fig.dpi)[1]/50)


#add grid
figs=(fig.get_size_inches()*fig.dpi).prod() #pixel area
h = g.get_yticks()
h=np.array([0]+h.tolist()+[max(h)+yticklabels]) #add zero line
w = g.get_xticks()
w=np.array([0]+w.tolist()+[max(w)+xticklabels]) #add zero line
g.hlines(h-0.5, w.max(),w.min(), linestyle='-', linewidth=0.5, color="black")
g.vlines(w-0.5,h.max(),h.min(), linestyle='-', linewidth=0.5, color="black")

plt.title("gene expression")
plt.ylabel("nutrient removal genes" )
plt.xlabel("genera" )
fig.savefig("gene_expression.png",dpi=400,bbox_inches="tight")
gdf.to_csv("gene expression.csv")

    #%%
#bias bar abosolute difference

fig = plt.figure()
defsize = fig.get_size_inches()*fig.dpi

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

patterns = [ "/" , "\\" , "|" , "-" , "+" , "x", "o", "O", ".", "*" ]

ind=pd.DataFrame([c.split(" ")[0] for c in gdf.columns]).drop_duplicates().set_index(0).index
dim=(len(ind),10)

av_allMPMG=allMPMG.loc[ind][::-1]/100#*100
av_allMP16S=allMP16S.loc[ind][::-1]/100
names=av_allMPMG.index.tolist()



floats=np.array(av_allMP16S)
colormap = cm.RdYlGn
normalize = mcolors.Normalize(vmin=np.min(floats), vmax=np.max(floats))
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)


sns.set_style("white")
fig, ax = plt.subplots(figsize=(dim[1]*0.5,dim[0]/3))
width=0.4
bottom=0
for x, i in enumerate(av_allMP16S):
    ax.barh(x, i, width, bottom, color=s_map.to_rgba(i), hatch=patterns[0],zorder=5)#, yerr=std_allMP16S[x] )
    ax.barh(x+0.4, av_allMPMG[x], width, bottom, color=s_map.to_rgba(av_allMPMG[x]), hatch=patterns[4],zorder=5)#, yerr=std_allMP16S[x] )
    
ax.invert_yaxis()
ax.set_yticks(np.arange(len(names))+0.25)
ax.set_yticklabels(names, Fontsize=9 )
ax.set_ylabel("genus", rotation=90)
ax.set_xlim([-0.2,0.1])
ax.set_ylim([-0.5,len(av_allMPMG)-0.2])



#hatches legend
circ1 = ax.barh(x, i, width, bottom, color='black',
                # alpha=0, 
                hatch=patterns[0])
circ2= ax.barh(x, i, width, bottom, color='black',
               # alpha=0, 
               hatch=patterns[4],zorder=0)
ax.legend([circ2,circ1],['MG-MP', '16S-MP'],loc=3)
ax.yaxis.tick_left()
ax.tick_params(axis='y',length=2)
plt.title(" Abundance Difference")
plt.xlabel("absolute abundance difference" )
plt.ylabel("genera" )
fig.savefig("genera_abundance_difference.png",dpi=400,bbox_inches="tight")
plt.show()

#%% bias bar percentage diference


#%% establish fold change order




c="/Volumes/Seagate_SSD/paper 3/third_paper_newest/Figures_main/4 Genus Metabolic /perc_change.xlsx"


cdf=pd.read_excel(c).fillna(0)
cdf=cdf[cdf["taxa"].str.startswith("g__")]
cdf=cdf.set_index("taxa")

ag=[i for i in allowed_genera if i in cdf.index]
cdf=cdf.loc[ag].fillna(0)

perc_MPMG=[]
perc_MP16S=[]

for plant in ["DXP","GW","SP"]:
    
    shMPMG=cdf[[plant+ " MP", plant+ " MG"]]
    shMP16S=cdf[[plant+ " MP", plant+ " S16"]]
    
    #shared
    shMPMG=shMPMG[(shMPMG>0).all(axis=1)]
    shMP16S=shMP16S[(shMP16S>0).all(axis=1)]
    
    #renormalized
    shMPMG=shMPMG/shMPMG.sum()*100
    shMP16S=shMP16S/shMP16S.sum()*100
    
    #percdiff
    perc_MPMG.append(shMPMG)
    perc_MP16S.append(shMP16S)
    
perc_MPMG=pd.concat(perc_MPMG,axis=1) 
MPMG_MP=perc_MPMG[[c for c in perc_MPMG.columns if "MP" in c]].sum(axis=1)/3
MPMG_MG=perc_MPMG[[c for c in perc_MPMG.columns if "MG" in c]].sum(axis=1)/3

perc_MP16S=pd.concat(perc_MP16S,axis=1)    
MP16S_MP=perc_MP16S[[c for c in perc_MP16S.columns if "MP" in c]].sum(axis=1)/3
MP16S_16S=perc_MP16S[[c for c in perc_MP16S.columns if "S16" in c]].sum(axis=1)/3
    #ok this doesnt work, you actually need to summ up the values for all the plants, not average the percentage change
    
#%% perc change
av_allMPMG=MPMG_MG-MPMG_MP
av_allMP16S=MP16S_16S-MP16S_MP

av_allMPMG=av_allMPMG.reset_index()
av_allMPMG.columns=["taxa","perc"]
av_allMP16S=av_allMP16S.reset_index()
av_allMP16S.columns=["taxa","perc"]

inddf=pd.DataFrame(ind.tolist(),columns=["taxa"])

av_allMPMG=inddf.merge(av_allMPMG,on="taxa",how="left").fillna(0).set_index("taxa")["perc"][::-1]
av_allMP16S=inddf.merge(av_allMP16S,on="taxa",how="left").fillna(0).set_index("taxa")["perc"][::-1]
names=av_allMPMG.index.tolist()



floats=np.array(av_allMP16S)
colormap = cm.RdYlGn
normalize = mcolors.Normalize(vmin=np.min(floats), vmax=np.max(floats))
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)


sns.set_style("white")
fig, ax = plt.subplots(figsize=(dim[1]*0.5,dim[0]/3))
width=0.4
bottom=0
for x, i in enumerate(av_allMP16S):
    #16S
    if i==2:
        plt.scatter(x,0.09,marker="o",s=40,c='r')
    elif i==-2:
        plt.scatter(x,0.09,marker="o",s=40,c='r')
    else:
        ax.barh(x, i, width, bottom, color=s_map.to_rgba(i), hatch=patterns[0],zorder=5)#, yerr=std_allMP16S[x] )
    
    #MG
    if  av_allMPMG[x]==2:
        plt.scatter(x,0.09,marker="o",s=40,c='r')
    elif av_allMPMG[x]==-2:
        plt.scatter(x,0.09,marker="o",s=40,c='r')
    
    ax.barh(x+0.4, av_allMPMG[x], width, bottom, color=s_map.to_rgba(av_allMPMG[x]), hatch=patterns[4],zorder=5)#, yerr=std_allMP16S[x] )


ax.invert_yaxis()
ax.set_yticks(np.arange(len(names))+0.5)
ax.set_yticklabels(names, Fontsize=9 )
ax.set_ylabel("genus", rotation=90)
ax.set_xlim([-20,15])
ax.set_ylim([-0.5,len(av_allMPMG)-0.2])



#hatches legend
circ1 = ax.barh(x, i, width, bottom, color='black',
                # alpha=0, 
                hatch=patterns[0])
circ2= ax.barh(x, i, width, bottom, color='black',
                # alpha=0, 
                hatch=patterns[4],zorder=0)
ax.legend([circ2,circ1],['MG-MP', '16S-MP'],loc=3)

plt.title(" Percentual Difference")
ax.yaxis.tick_left()
ax.tick_params(axis='y',length=2)
plt.xlabel("absolute abundance difference" )
plt.ylabel("genera" )
fig.savefig("genera_percentual_difference.png",dpi=400,bbox_inches="tight")
plt.show()

#%% percentual
av_allMPMG=(MPMG_MG-MPMG_MP).divide(((MPMG_MG+MPMG_MP)/2))
av_allMP16S=(MP16S_16S-MP16S_MP).divide(((MP16S_16S+MP16S_MP)/2))

av_allMPMG=av_allMPMG.reset_index()
av_allMPMG.columns=["taxa","perc"]
av_allMP16S=av_allMP16S.reset_index()
av_allMP16S.columns=["taxa","perc"]

inddf=pd.DataFrame(ind.tolist(),columns=["taxa"])

av_allMPMG=inddf.merge(av_allMPMG,on="taxa",how="left").fillna(0).set_index("taxa")["perc"][::-1]
av_allMP16S=inddf.merge(av_allMP16S,on="taxa",how="left").fillna(0).set_index("taxa")["perc"][::-1]
names=av_allMPMG.index.tolist()



floats=np.array(av_allMP16S)
colormap = cm.RdYlGn
normalize = mcolors.Normalize(vmin=np.min(floats), vmax=np.max(floats))
s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)


sns.set_style("white")
fig, ax = plt.subplots(figsize=(dim[1]*0.5,dim[0]/3))
width=0.4
bottom=0
for x, i in enumerate(av_allMP16S):
    
    #16S
    if i==2:
        plt.scatter(x,0.09,marker="o",s=40,c='r')
    elif i==-2:
        plt.scatter(x,0.09,marker="o",s=40,c='r')
    else:
        ax.barh(x, i, width, bottom, color=s_map.to_rgba(i), hatch=patterns[0],zorder=5)#, yerr=std_allMP16S[x] )
    
    #MG
    if  av_allMPMG[x]==2:
        plt.scatter(x,0.09,marker="o",s=40,c='r')
    elif av_allMPMG[x]==-2:
        plt.scatter(x,0.09,marker="o",s=40,c='r')
    
    ax.barh(x+0.4, av_allMPMG[x], width, bottom, color=s_map.to_rgba(av_allMPMG[x]), hatch=patterns[4],zorder=5)#, yerr=std_allMP16S[x] )


    
ax.invert_yaxis()
ax.set_yticks(np.arange(len(names))+0.5)
ax.set_yticklabels(names, Fontsize=9 )
ax.set_ylabel("genus", rotation=90)
ax.set_xlim([-2,2])
ax.set_ylim([-0.5,len(av_allMPMG)-0.2])



#hatches legend
circ1 = ax.barh(x, i, width, bottom, color='black',
                # alpha=0, 
                hatch=patterns[0])
circ2= ax.barh(x, i, width, bottom, color='black',
                # alpha=0, 
                hatch=patterns[4],zorder=0)
ax.legend([circ2,circ1],['MG-MP', '16S-MP'],loc=3)

plt.title(" Percentual Difference")
ax.yaxis.tick_left()
ax.tick_params(axis='y',length=2)
plt.xlabel("percentual abundance difference" )
plt.ylabel("genera" )
fig.savefig("genera_percentual_difference.png",dpi=400,bbox_inches="tight")
plt.show()

