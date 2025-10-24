# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 12:45:06 2025

@author: e_kle
"""
#https://stackoverflow.com/questions/17543359/drawing-lines-between-two-plots-in-matplotlib
import numpy as np
import pandas as pd
import os
from pathlib import Path

files=[
"E:/Data/Data_analysis/Skeletons_OrbiSIMS/MFP_GOld_Output/gold_reference_LMIG_orbi_(-)1_all_peaks.csv",
"E:/Data/Data_analysis/Skeletons_OrbiSIMS/MFP_GOld_Output/gold_reference_GCIB_orbi_(-)1_all_peaks.csv"]

cpre=os.path.commonprefix(files)
csuf=os.path.commonprefix([i[::-1] for i in files])[::-1]


import matplotlib.pyplot as plt
import matplotlib


rows,cols=len(files),1
fig,axes=plt.subplots(rows,cols,figsize=(10,rows*4))
if (rows*cols)==1: axes=np.array([axes])
axes = axes.flatten()

test=[]

ds=[]
for i,f in enumerate(files): 
    d=pd.read_csv(f)[["Mass","Apex"]]
    d["ax"]=i
    ds.append(d)
    
    x,y=d["Mass"],d["Apex"]
    axes[i].scatter(x,y,color="black")    
    test.append([x[0],y[0]])
    axes[i].set_xlim(0,1000)
    #axes[i].set_ylabel(Path(f).stem.replace('_all_peaks',""))
    
    axes[i].set_ylabel(f.replace(cpre,"").replace(csuf,""))
    
    #find shared parts between stems
    

    
    #clean frames and ticks
    if i<(len(files)-1):
          axes[i].spines['bottom'].set_visible(False)
          axes[i].get_xaxis().set_ticklabels([])
    if i: axes[i].spines['top'].set_visible(False)


plt.suptitle(Path(cpre).stem.strip("_")+csuf.replace("_all_peaks.csv",""))

ds=pd.concat(ds).sort_values(by="Mass").reset_index(drop=True)

#find connections
#each point can only be connected once between files
#the connection should be between the most intense points
#kde against itself
#generate format: 
    #ix, mass, Apex, ax, 
    #connecting ix, mass, Apex, ax
    #unique ix,ix_c combination
    #total Apex
    
#groupby ix,ix_c compbination pick highest Apex
from sklearn.neighbors import KDTree
ppm=5
m=ds.Mass.values
tree = KDTree(m.reshape(1,-1).T, leaf_size=20) 
t=tree.query_radius(m.reshape(1,-1).T,r=m*5/1e6)
tix=np.vstack([np.hstack([np.repeat(ix,len(i)) for ix,i in enumerate(t)]),np.hstack(t)]).T

#trim unique combinations 
tix=tix[tix[:,0]!=tix[:,1]]
tix=tix[ds.ax.loc[tix[:,0]].values!=ds.ax.loc[tix[:,1]].values]

#ac=np.unique(ds.loc[np.unique(tix)].ax)
tdf=pd.DataFrame(tix,columns=["a","b"])
tdf["s"]=[str(set(i)) for i in tdf[["a","b"]].values]
tdf["t"]=ds.loc[tdf.a,"Apex"].values+ds.loc[tdf.b,"Apex"].values
best=tdf.sort_values(by=["s","t"],ascending=False).groupby("s",sort=False).nth(0) #most intense connections

#plot connections
fig.canvas.draw()
transFigure = fig.transFigure.inverted()
for r in best[["a","b"]].values:
    a,b=ds.loc[r[0]],ds.loc[r[1]]
    c1 = transFigure.transform(axes[int(a.ax)].transData.transform(a.values[:2]))
    c2 = transFigure.transform(axes[int(b.ax)].transData.transform(b.values[:2]))
    fig.lines.append( matplotlib.lines.Line2D((c1[0],c2[0]),(c1[1],c2[1]),
                                    transform=fig.transFigure,linestyle="--",color="grey",linewidth=1))

#%%
#add labels to top x most abundant pairs
#add y-axis labels

#add MFP?
#trim isotopes?


