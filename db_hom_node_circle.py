# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 13:14:19 2021

@author: hbckleikamp
"""
#%% clear variables

# try:
#     from IPython import get_ipython
#     get_ipython().magic('clear')
#     get_ipython().magic('reset -f')
# except:
#     pass

#%% change directory to script directory
import os
from pathlib import Path
import sys
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())

#%% standard Modules
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt



from itertools import chain, groupby
from collections import Counter
from openpyxl import load_workbook 

import Bio
from Bio import SeqIO
import seaborn as sns; sns.set()
import time

from collections import OrderedDict

#%% Read GTDB ssu file
ranks=["superkingdom","phylum","class","order","family","genus","species"]
prefix=["d__","p__","c__","o__","f__","g__","s__"]
file="/Volumes/Seagate_SSD/third_paper new/Annotation/16S/2_Prepping_GTDB_16S/ssu_all.tar"
records=pd.DataFrame([[str(i.id),str(i.seq)] for i in Bio.SeqIO.parse(file,"fasta")],columns=["Accession","seqs"])
records["Accession"]=records["Accession"].apply(lambda x: x.split("~")[0])
records["length"]=records["seqs"].apply(len)

#%%

# #records.hist("length")
# plt.figure()
# plt.hist(records["length"].tolist(),bins=100)
# plt.xlim([0,1700])

#%% merge with right on metadata GTDB taxonomies:
    # -not in 16S
    # -fl
    # -nfl

taxfiles=["/Volumes/Seagate_SSD/third_paper new/Data resources/Public databases/GTDB/gtdb metadata/bac120_taxonomy.tsv",
"/Volumes/Seagate_SSD/third_paper new/Data resources/Public databases/GTDB/gtdb metadata/ar122_taxonomy.tsv"]

tdf=[]
for t in taxfiles:
    tdf.append(pd.read_csv(t,header=None,sep="\t"))#["Accession"]+ranks)
tdf=pd.concat(tdf)
tdf.columns=["Accession","lineage"]
tdf[ranks]=tdf["lineage"].str.rsplit(";",expand=True)
tdf=tdf.merge(records,how="left",on="Accession").fillna(0)


not_in_16S=tdf[tdf["length"]==0].groupby("lineage").size().reset_index()
not_in_16S.columns=["lineage","Count"]
not_in_16S[ranks]=not_in_16S["lineage"].str.rsplit(";",expand=True)
not_in_16S=not_in_16S[ranks+["Count"]]

not_full_length=tdf[(tdf["length"]>0) & (tdf["length"]<1200)].groupby("lineage").size().reset_index()
not_full_length.columns=["lineage","Count"]
not_full_length[ranks]=not_full_length["lineage"].str.rsplit(";",expand=True)
not_full_length=not_full_length[ranks+["Count"]]

full_length=tdf[tdf["length"]>=1200].groupby("lineage").size().reset_index()
full_length.columns=["lineage","Count"]
full_length[ranks]=full_length["lineage"].str.rsplit(";",expand=True)
full_length=full_length[ranks+["Count"]]





#%% make phylum colored trees?

class Node:
    def __init__(self,name,rank,count,xcor,ycor):
        self.name=name
        self.rank=rank
        self.count=count
        self.xcor=xcor
        self.ycor=ycor
        #not each node has a parent
        


#%%get all phylum names for each rank to make a shared color scheme



    

top8=[i[0] for i in Counter(tdf["phylum"].tolist()).most_common(20)]
colordict={}
cmap = sns.color_palette("tab20")
[colordict.update({i:cmap[c]}) for c,i in enumerate(top8)]



#%% nested list
ranks=["root","superkingdom","phylum","class","order","family","genus","species"]

for df in [not_in_16S,not_full_length,full_length]:
    break
df=df.iloc[0:1000]
df["root"]=["root"]*len(df)


#fill empty nodes
emptyinds=np.where(df.to_numpy()=="")
for c,i in enumerate(emptyinds[0]):
    df.iloc[emptyinds[0][c],emptyinds[1][c]]="empty_"+str(c)




#create nodes

nodes=[]
node_names=[]
rank_groups=[]

space_top=10

for r_ix,rank in enumerate(ranks):

    counts=df.groupby(rank)["Count"].sum()/df["Count"].sum()*100
    names=df[rank].drop_duplicates()
    rank_groups.append(names)
    ycors=np.linspace(0,space_top,len(names))

    

    
    
    for n_ix,name in enumerate(names):
        
        print(n_ix)
        if name!="":
            count=counts.loc[name]
            if "empty_" in name: count=0
            
            #define node
            mynode=Node(name,rank,count,
                            xcor=r_ix,
                            ycor=ycors[n_ix])
                
            #add parents
            if r_ix: 
                mynode.parent=df.loc[df[ranks[r_ix]]==name,ranks[r_ix-1]].values[0]
                          

            #add children
            if r_ix<len(ranks)-1: 
                children=np.unique(df.loc[df[ranks[r_ix]]==name,ranks[r_ix+1]].values).tolist()
                mynode.children=[i for i in children if i !=""]       
            
        
            #add color based on phylum rank
            mynode.color='grey'
            rname=df.loc[df[rank]==name,"phylum"].tolist()[0]
            if rname in top8:                        
                mynode.color=colordict.get(rname)
                
            
            node_names.append(mynode.name)
            nodes.append(mynode)
                    
                    

            
node_dict =  OrderedDict(zip(node_names, nodes))
#test later if using a dict or ordereddict changes output?



#center parents to children
for ix,r in enumerate(rank_groups[::-1]):
    
    if ix:
        for i in r:
            node=node_dict.get(i)
            ycors=[node_dict.get(c).ycor for c in node.children]
            node.ycor=(max(ycors)+min(ycors))/2
            



#project on circle circumference
cors=[]
for k,node in node_dict.items():
    angle=360*node.ycor/space_top
    radius=node.xcor
    node.ycor=radius*math.sin(math.radians(angle))
    node.xcor=radius*math.cos(math.radians(angle))

#define parent cors
for ix,r in enumerate(rank_groups):
      if ix:
         for i in r:
            node=node_dict.get(i)
            node.parentycor=node_dict.get(node.parent).ycor
            node.parentxcor=node_dict.get(node.parent).xcor
    




#%%   

#scale figsize (ylen) with len of vals?
figsize=(60,40)  #xlen, ylen
plt.figure(figsize=figsize) #scale thickness with figsize
for node in nodes:
 
    
    if hasattr(node, 'parent'):
        
        #plot hanging line
        # point1=[node.xcor,node.ycor]#child
        # point2=[node.xcor-1,node.parentycor] #parent
        # a = (point2[1] - point1[1])/(np.cosh(point2[0]) - np.cosh(point1[0]))
        # b = point1[1] - a*np.cosh(point1[0])
        # x = np.linspace(point1[0], point2[0], 100)
        # y = a*np.cosh(x) + b
        
        #plot straight line
        x=[node.xcor,node.parentxcor]
        y=[node.ycor,node.parentycor]        
        plt.plot(x, y, c=node.color
                  ,linewidth=node.count*0.4*(figsize[1]/5), alpha=1)

        

for node in nodes:
    
    if node.rank!="superkingdom_name":
        plt.plot(node.xcor,node.ycor,".",color=node.color
                 ,markersize=node.count*figsize[1]/5, alpha=1)

    
  
#plt.xlim([-0.2,7.2])

plt.axis('off')
plt.grid(b=None)
 
plt.show()




        
#%% legend plot
    
    

cdf=pd.DataFrame(colordict.items(),columns=["name","color"]).set_index("name")


#cdf=cdf.T[podf].T

plt.figure()
markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in cdf["color"].tolist()]
plt.axis('off')
plt.grid(b=None)

plt.legend(markers, [i[3:] for i in cdf.index.tolist()], numpoints=1,facecolor='white')
plt.show()
