# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:03:30 2024

@author: hkleikamp
"""

import pandas as pd

#%% files

kf="E:/Novolign/0508 Revisions/2024-08-05keg.tsv"


kos=["E:/Novolign/0508 Revisions/KOs_DB_searching_proteins/KOs_DB_searching_proteins/user_ko_rep01.txt",
"E:/Novolign/0508 Revisions/KOs_DB_searching_proteins/KOs_DB_searching_proteins/user_ko_rep02.txt"]

psms=[
"E:/Novolign/0508 Revisions/PEAKS_DB_SEARCHING_OUTPUT_UProt_DB/PEAKS_DB_SEARCHING_OUTPUT_UProt_DB/MP_DN_Study_30072024_DDA01_UP man DB_15_AST_cea_30min_DDA01_Rep01/DB search psm.csv",
"E:/Novolign/0508 Revisions/PEAKS_DB_SEARCHING_OUTPUT_UProt_DB/PEAKS_DB_SEARCHING_OUTPUT_UProt_DB/MP_DN_Study_30072024_DDA01_UP man DB_15_AST_cea_30min_DDA01_Rep02/DB search psm.csv"]


#%%



kdf=pd.read_csv(kf,sep="\t")

#clean up for figure
k=kdf[['cat1', 'cat2', 'cat3',"KO"]]
k["cat3"]=k["cat3"].str.split("[").apply(lambda x: x[0])

for i in ['cat1', 'cat2', 'cat3']:
    k[i]=k[i].apply(lambda x: x[6:]).str.strip()


cats=['Metabolism', 'Genetic Information Processing',  'Environmental Information Processing', 'Cellular Processes']
k=k[k['cat1'].isin(cats)]


#%%




for ix,i in enumerate(kos):




    with open(i) as f:
        lines=f.readlines()
    
    fdf=pd.DataFrame(lines,columns=["lines"])
    fdf[["acc","KO"]]=fdf["lines"].str.replace("\n","").str.rsplit("\t",expand=True).fillna("")
    
    pdf=pd.read_csv(psms[ix])
    
    
    #Im not doing any FDR filtering or anything just going straight to the accessions
    adf=pdf["Accession"].str.split(":").explode().fillna("").to_frame("acc")
    s=adf.groupby("acc").size().to_frame("count").reset_index()
    
    d=fdf[["acc","KO"]].merge(s,on="acc",how="right").merge(k,how="inner",on="KO")
    
    
    groups=d.groupby("cat1")
    
    #%% Plotly method
    
    import plotly.express as px
    import plotly.io as pio
    import plotly.graph_objects as go
    pio.renderers.default='browser'
    
    for n,g in groups:
        fig = px.sunburst(g,path=['cat2', 'cat3'], values='count')
        fig.show()
    #%%
    
    #stacked pie method (fake doughnut)
    #https://stackoverflow.com/questions/44153457/double-donut-chart-in-matplotlib

 
    
    # al='abcdefghijklmnopqrstuvwxyz'
    
    # import matplotlib.pyplot as plt
    # import numpy as np
    
    # fig, ax = plt.subplots()
    # ax.axis('equal')
    # width = 0.3
    
    # # Outer ring
    
    # r1=g.groupby("cat2")["count"].sum()
    # l1=[al[i].upper() for i in range(len(r1))]
    
    # r2=g.groupby(["cat2","cat3"])["count"].sum()
    
    # dc=dict(list(zip(r1.index.tolist(),l1)))
    
    # let=[dc.get(i) for i in r2.reset_index()["cat2"]]
    
    # l2=(pd.Series(let)+np.hstack([np.arange(len(i2)) for i1,i2 in r2.groupby(["cat2"])]).astype(str)).tolist()
       
        
    
    
    # cm = plt.get_cmap("tab20c")
    # cout = cm(np.arange(3)*4)
    # pie, _ = ax.pie(r2.values, radius=1, labels=l2, colors=cout)
    # plt.setp(pie, width=width, edgecolor='white')
    
 
    
    # # Inner ring
    # cin = cm(np.array([1,2,5,6,9,10]))
    # labels = list(map("".join, zip(list("aabbcc"),map(str, [1,2]*3))))
    # pie2, _ = ax.pie(r1.values, radius=1-width, labels=l1,
    #                   labeldistance=1.2, colors=cin)
    # plt.setp(pie2, width=width, edgecolor='white')
    # plt.show()

    