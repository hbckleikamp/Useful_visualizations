# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 21:29:15 2022

@author: ZR48SA
"""

import pandas as pd
from pathlib import Path

# define lineage ranks
ranks=["superkingdom","phylum","class","order","family","genus","species"] 
gtdb_ranks=["GTDB_"+r for r in ranks]
rank_names=[r+"_name" for r in ranks]
rank_ids=[r+"_id" for r in ranks]


# files=[
# "C:/Comet/target_C6_1_JSP_JLN_UPLIFT_annotated.csv",
# "C:/Comet/target_C6_2_JSP_JLN_UPLIFT_annotated.csv",
# "C:/Comet/target_C24_1_JSP_JLN_UPLIFT_annotated.csv",
# "C:/Comet/target_C24_2_JSP_JLN_UPLIFT_annotated.csv",
# "C:/Comet/target_Re1_JSP_JLN_UPLIFT_annotated.csv",
# "C:/Comet/target_Ox2_JSP_JLN_UPLIFT_annotated.csv",
# "C:/Comet/target_Re2_JSP_JLN_UPLIFT_annotated.csv",
# "C:/Comet/target_Ox1_JSP_JLN_UPLIFT_annotated.csv"

# ]


files=["D:/Analysis_Paola2/merged/merged_target_T0_1_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T0_2_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T0_A1_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T0_A2_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T1_1_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T1_2_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T1_A1_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T1_A2_HK_JLN_UPLIFT.tsv"]



dfs=[]
for f in files:
    df=pd.read_csv(f)
    #dfs.append(df[["spectrum_reference"]+ranks].drop_duplicates().groupby("genus").size())
    dfs.append(df[["First Scan"]+gtdb_ranks].drop_duplicates().groupby("GTDB_genus").size())
dfs=pd.concat(dfs,axis=1)
dfs.columns=[Path(file).stem for file in files]


dfs["sum"]=dfs.sum(axis=1)

dfs=dfs.sort_values(by="sum")
dfs=dfs.iloc[-30:]
dfs=dfs/dfs.sum()
dfs=dfs.round(2)
dfs=dfs[::-1]


rallcomb=dfs[[Path(file).stem for file in files]]
filenames=["T0_1",
            "T0_2",
            "T0_A1",
            "T0_A2",
            "T1_1",
            "T1_2",
            "T1_A1",
            "T1_A2"]

rallcomb.columns=filenames



import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("white")
fig, ax = plt.subplots(figsize=(4, 12))
# yticklabels=3
# xticklabels=3
g=sns.heatmap(rallcomb,
               cmap="YlGnBu",
               #center=20,
               #cmap="vlag",
               mask=(rallcomb==0),
               # yticklabels=yticklabels,
               robust=True,
               #vmin=0, vmax=40,
               cbar='right'
               )




g.set_yticklabels([str(i).split("'")[1].split(" ")[0].split("__")[1] for i in g.get_yticklabels()],
                  fontsize=(fig.get_size_inches()*fig.dpi)[0]/19,
                  va='top')
g.set_xticklabels(g.get_xticklabels(), 
                  rotation=90, 
                  fontsize = (fig.get_size_inches()*fig.dpi)[1]/50)

fig.savefig('hm.png',dpi=600,bbox_inches="tight")

#%%
files=["D:/Analysis_Paola2/calisp_output/T0_1_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T0_2_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T0_A1_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T0_A2_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T1_1_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T1_2_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T1_A1_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T1_A2_HK_JLN_UPLIFT.feather"]



cdfs=[]
for ix,file in enumerate(files):
    df=pd.read_feather(file)
    #df.ratio_na.plot.density()

    #df=pd.read_feather(files[6])
    df=df[df["bins"]!="decoy"]
    cdf=pd.DataFrame(df["ratio_na"]) #.plot.hist(bins=40)
    cdf=cdf[cdf<cdf.quantile(0.95)*1.4]
    
    print(cdf.median())
    
    
    cdf["g"]=[filenames[ix]]*len(cdf)
    
    cdf=cdf
    cdfs.append(cdf)

# #boxplot

# dat=pd.concat(cdfs)
# import seaborn as sns
# sns.set_style("white")
# fig,ax=plt.subplots()

# #ax=rdf[["role_genes","taxon_rank"]].boxplot(by='role_genes')

# ax = sns.boxplot(x="g", y="ratio_na",data=dat)


#used sns ridgeplot: https://seaborn.pydata.org/examples/kde_ridgeplot.html


#seaborn components used: set_theme(), cubehelix_palette(), FacetGrid

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})


df=pd.concat(cdfs).reset_index(drop=True)
df.columns=['x','g']




# Initialize the FacetGrid object
pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
pal = [i+[0.5] for i in sns.cubehelix_palette(10, rot=-.25, light=.6)]

fig,ax=plt.subplots()
g = sns.FacetGrid(df, row="g", hue="g", aspect=15, height=.5, palette=pal)

# Draw the densities in a few steps
g.map(sns.kdeplot, "x",
      bw_adjust=.5, clip_on=False,
      fill=True, alpha=1, linewidth=1.5)
g.map(sns.kdeplot, "x", clip_on=False, color="w", lw=2, bw_adjust=.5)

# passing color=None to refline() uses the hue mapping
g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)


# Define and use a simple function to label the plot in axes coordinates
def label(x, color, label):
    ax = plt.gca()
    ax.text(-.1, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)


g.map(label, "x")

# Set the subplots to overlap
g.figure.subplots_adjust(hspace=-.25)

# Remove axes details that don't play well with overlap
g.set_titles("")
g.set(yticks=[], ylabel="")
g.set(xlabel="labelling ratio")
g.despine(bottom=True, left=True)
g.set(xlim=(0, 0.03))

#fig = ax.get_figure()
g.savefig('na.png',dpi=600)


#%%

import pandas as pd
from pathlib import Path

# define lineage ranks
ranks=["superkingdom","phylum","class","order","family","genus","species"] 
rank_names=[r+"_name" for r in ranks]
rank_ids=[r+"_id" for r in ranks]


files=[
"D:/Analysis_Paola2/merged/merged_target_T0_1_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T0_2_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T0_A1_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T0_A2_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T1_1_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T1_2_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T1_A1_HK_JLN_UPLIFT.tsv",
"D:/Analysis_Paola2/merged/merged_target_T1_A2_HK_JLN_UPLIFT.tsv"]

cps=["D:/Analysis_Paola2/calisp_output/T0_1_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T0_2_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T0_A1_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T0_A2_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T1_1_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T1_2_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T1_A1_HK_JLN_UPLIFT.feather",
"D:/Analysis_Paola2/calisp_output/T1_A2_HK_JLN_UPLIFT.feather"]

fig,ax=plt.subplots()
ecs=[]
gos=[]
genera=[]
iprs=[]
for ix,file in enumerate(files):




    fdf=pd.read_csv(file) 
    cdf=pd.read_feather(cps[ix]) 
    cdf["peptide"]=cdf["peptide"].apply(lambda x: x.split()[0])
    
    fdf=fdf.merge(cdf,left_on="Sequence",right_on="peptide",how="left")
    
    
    ecs.append(fdf[["sequence","ratio_fft","ec"]].drop_duplicates().groupby("ec")["ratio_fft"].mean())
    gos.append(fdf[["sequence","ratio_fft","go"]].drop_duplicates().groupby("go")["ratio_fft"].mean())
    iprs.append(fdf[["sequence","ratio_fft","ipr"]].drop_duplicates().groupby("ipr")["ratio_fft"].mean())
    genera.append(fdf[["sequence","ratio_fft","genus"]].drop_duplicates().groupby("genus")["ratio_fft"].mean())


ecs=pd.concat(ecs,axis=1)
ecs.columns=filenames
ecs["sum"]=ecs.sum(axis=1)
ecs=ecs.sort_values(by="sum")
ecs=ecs.iloc[-30:]
ecs=ecs[::-1]
ecs=ecs[filenames]

sns.set_style("white")
fig, ax = plt.subplots(figsize=(4, 12))
# yticklabels=3
# xticklabels=3
g=sns.heatmap(ecs,
               cmap="YlGnBu",
               #center=20,
               #cmap="vlag",
               mask=(ecs==0),
               # yticklabels=yticklabels,
               robust=True,
               #vmin=0, vmax=40,
               cbar='right'
               )




# g.set_yticklabels([str(i).split("'")[1].split(" ")[0].split("__")[1] for i in g.get_yticklabels()],
#                   fontsize=(fig.get_size_inches()*fig.dpi)[0]/19,
#                   va='top')
g.set_xticklabels(g.get_xticklabels(), 
                  rotation=90, 
                  fontsize = (fig.get_size_inches()*fig.dpi)[1]/50)

fig.savefig('echm.png',dpi=600,bbox_inches="tight")


gos=pd.concat(gos,axis=1)
gos.columns=filenames
gos["sum"]=gos.sum(axis=1)
gos=gos.sort_values(by="sum")
gos=gos.iloc[-30:]
gos=gos[::-1]
gos=gos[filenames]

sns.set_style("white")
fig, ax = plt.subplots(figsize=(4, 12))
# yticklabels=3
# xticklabels=3
g=sns.heatmap(gos,
               cmap="YlGnBu",
               #center=20,
               #cmap="vlag",
               mask=(gos==0),
               # yticklabels=yticklabels,
               robust=True,
               #vmin=0, vmax=40,
               cbar='right'
               )




# g.set_yticklabels([str(i).split("'")[1].split(" ")[0].split("__")[1] for i in g.get_yticklabels()],
#                   fontsize=(fig.get_size_inches()*fig.dpi)[0]/19,
#                   va='top')
g.set_xticklabels(g.get_xticklabels(), 
                  rotation=90, 
                  fontsize = (fig.get_size_inches()*fig.dpi)[1]/50)

fig.savefig('gohm.png',dpi=600,bbox_inches="tight")


iprs=pd.concat(iprs,axis=1)
iprs.columns=filenames
iprs["sum"]=iprs.sum(axis=1)
iprs=iprs.sort_values(by="sum")
iprs=iprs.iloc[-30:]
iprs=iprs[::-1]
iprs=iprs[filenames]

sns.set_style("white")
fig, ax = plt.subplots(figsize=(4, 12))
# yticklabels=3
# xticklabels=3
g=sns.heatmap(iprs,
               cmap="YlGnBu",
               #center=20,
               #cmap="vlag",
               mask=(iprs==0),
               # yticklabels=yticklabels,
               robust=True,
               #vmin=0, vmax=40,
               cbar='right'
               )




# g.set_yticklabels([str(i).split("'")[1].split(" ")[0].split("__")[1] for i in g.get_yticklabels()],
#                   fontsize=(fig.get_size_inches()*fig.dpi)[0]/19,
#                   va='top')
g.set_xticklabels(g.get_xticklabels(), 
                  rotation=90, 
                  fontsize = (fig.get_size_inches()*fig.dpi)[1]/50)

fig.savefig('iprhm.png',dpi=600,bbox_inches="tight")


genera=pd.concat(genera,axis=1)
genera.columns=filenames
genera["sum"]=genera.sum(axis=1)
genera=genera.sort_values(by="sum")
genera=genera.iloc[-30:]
genera=genera[::-1]
genera=genera[filenames]

sns.set_style("white")
fig, ax = plt.subplots(figsize=(4, 12))
# yticklabels=3
# xticklabels=3
g=sns.heatmap(genera,
               cmap="YlGnBu",
               #center=20,
               #cmap="vlag",
               mask=(genera==0),
               # yticklabels=yticklabels,
               robust=True,
               #vmin=0, vmax=40,
               cbar='right'
               )




# g.set_yticklabels([str(i).split("'")[1].split(" ")[0].split("__")[1] for i in g.get_yticklabels()],
#                   fontsize=(fig.get_size_inches()*fig.dpi)[0]/19,
#                   va='top')
g.set_xticklabels(g.get_xticklabels(), 
                  rotation=90, 
                  fontsize = (fig.get_size_inches()*fig.dpi)[1]/50)

fig.savefig('generahm.png',dpi=600,bbox_inches="tight")


#%%
    
# sequence)
# ratio_fft
# 	0
# 35	ec
# 36	go
# 37	ipr

#%%

dfs=[]
for ix,f in files:
    df=pd.read_csv(f)
    dfs.append(df[["spectrum_reference"]+ranks].drop_duplicates().groupby("genus").size())
dfs=pd.concat(dfs,axis=1)
dfs.columns=[Path(file).stem for file in files]
