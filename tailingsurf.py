#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 20:33:50 2021

@author: hugokleikamp
"""
#%%
import sys
from inspect import getsourcefile
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend import Legend

import random, re, requests
import threading, time, string
import pickle, math

from pathlib import Path #path package is better?
from path import Path as path

from scipy.optimize import curve_fit
from itertools import chain, groupby #look at module inconsistency
import itertools
from collections import Counter

import datetime
import dateutil.parser as dparser
import ftputil, urllib, gzip, zipfile, shutil
from openpyxl import load_workbook #for kronaplots
# import Bio
# from Bio import SeqIO                         #for homology search
import proteoclade, sqlite3, hashlib #for proteoclade annotation

#sunburst plot   
import plotly
import plotly.express as px
import plotly.io as pio
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





taxdf=pd.read_excel("/Volumes/Seagate_SSD/third_paper/outputs for figures/rank_quant/MG_quant_DXP.xlsx")
countcols=[i for i in taxdf.columns if "count" in i]


#%%



#%%
        
#plot sploping
#sunburst plot   
import plotly
import plotly.express as px
import plotly.io as pio

import plotly.graph_objects as go

import numpy as np
pio.renderers.default='browser'
#pio.renderers.default='svg'
#
#%%




cmap = plt.get_cmap("tab10")
custom_colorscale=[[0.0, 'rgb' + str(cmap(2)[0:3])],
                  [0.5, 'rgb' + str(cmap(0)[0:3])],
                  [1, 'rgb' + str(cmap(4)[0:3])]]





ranks=["superkingdom","phylum","class","order","family","genus","species"][1:]
       
      
folder="/Volumes/Seagate_SSD/third_paper new/Annotation/Output files/rank_quant/"    
files=os.listdir(folder)
       
for sample in ["DXP","GW","SP"]:
    
    
    samplefiles=[f for f in files if sample in f]
    MPfile=[f for f in samplefiles if "MPDB" in f][0]
    MGfile=[f for f in samplefiles if "MG" in f][0]
    S16file=[f for f in samplefiles if "S16" in f][0]
    MP=pd.read_excel(str(Path(folder,MPfile)))[countcols].values[:,1:]
    MG=pd.read_excel(str(Path(folder,MGfile)))[countcols].values[:,1:]
    S16=pd.read_excel(str(Path(folder,S16file)))[countcols].values[:,1:]
    
    
    
    fig=go.Figure(data=[go.Surface(z=MP,                           
                                surfacecolor=np.ones(MP.shape)*0,
                                opacity=.2, 
                                cmin=0,
                                cmax=1,
                                colorscale=custom_colorscale,
                                showscale=False
                                ),
                  
                         go.Surface(z=MG,                           
                                surfacecolor=np.ones(MG.shape)*0.5,
                                opacity=.2, 
                                cmin=0,
                                cmax=1,
                                colorscale=custom_colorscale,
                                showscale=False
                                ),
                         
                         
                        go.Surface(z=S16,                           
                                surfacecolor=np.ones(S16.shape)*1,
                                opacity=.2, 
                                cmin=0,
                                cmax=1,
                                colorscale=custom_colorscale,
                                showscale=False
                                )
                       
                          
                          
                          ]                            
                  )
    
    
    
    name = sample
    
    camera = dict(
        up=    dict(x=0, y=0,   z=1),
        center=dict(x=0, y=0,   z=0),
        eye=   dict(x=2, y=1.2, z=0.3)
    )
    
    
    
    fig.update_layout(scene_camera=camera, title=name)  #, yaxis=yaxis)
    
    fig.update_layout(
        scene = dict(
            xaxis_title='ranks',
            yaxis_title='taxa',
            zaxis_title='percent abundance',
            
            zaxis = dict(range=[0,75]),
            yaxis = dict(range=[0,130]),
    
            
            xaxis = dict(
            tickmode = 'array',
            tickvals = np.arange(len(ranks)),
            ticktext = ranks
            
    
            
            )))
    
    fig.show()


