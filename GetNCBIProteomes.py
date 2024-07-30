# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:52:48 2024

@author: hkleikamp
"""

#%%

import pandas as pd
import numpy as np
from io import StringIO
import os
from pathlib import Path

from urllib.request import urlretrieve
import gzip, zipfile, shutil, tarfile, requests, time
#%%
def extract(path): #path to folder

    #recursive extraction            
    while any([f.endswith((".zip",".gz",".tar")) for f in os.listdir(path)]):
        
        for f in os.listdir(path):
            
            i=str(Path(path,Path(f)))
            o=str(Path(path,Path(f).stem))
            
            if f[0].isalnum() and f.endswith(".zip"):
                print("extracting "+f)
                with zipfile.ZipFile(i, 'r') as zip_ref:
                    zip_ref.extractall(path)
                if os.path.exists(i): os.remove(i)

            if f[0].isalnum() and f.endswith(".gz"):
                print("extracting "+f)
                with gzip.open(i,'rb') as f_in:
                    with open(o,'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                if os.path.exists(i): os.remove(i)
                
            if f[0].isalnum() and f.endswith(".tar"):
                print("extracting "+f)
                tar = tarfile.open(i, "r:")
                tar.extractall(path)
                tar.close()
                if os.path.exists(i): os.remove(i)
    


#%%


#read input file
file="C:/Users/hkleikamp/Downloads/SI3_List of organisms selected for building DB48 database.xlsx"
df=pd.read_excel(file,engine="openpyxl")
c=df.iloc[1,:]
df=df.iloc[2:]
df.columns=c
taxids=df["Species Taxid"].values.astype(str)

#Get proteomes

url_l="https://rest.uniprot.org/proteomes/stream?compressed=false&fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd%2Cgenome_assembly&format=tsv&query=%28organism_id%3A"
url_r="%29"



proteomes=[]
for i in taxids:
    
    url=url_l+i+url_r
    print(i)
    print(url)
    r=requests.get(url)
    proteomes.append(pd.read_csv(StringIO(r.text),sep="\t"))
  
proteomes=pd.concat(proteomes).reset_index(drop=True)

#parse busco
proteomes["BUSCO"]=proteomes["BUSCO"].str.replace("[",",",regex=False)
proteomes["BUSCO"]=proteomes["BUSCO"].str.replace("]","",regex=False)
busco=(proteomes["BUSCO"]+["%"]).str.rsplit(",",expand=True).applymap(lambda x: x[2:-1]).astype(float)
busco.columns=["completeness","single","duplicated","fragmented","missing","genes"]

#pick most complete proteome
df=pd.concat([proteomes,busco],axis=1)
df=df.sort_values(by=["Organism Id","completeness"],ascending=False)
bdf=df.groupby("Organism Id",sort=False).nth(0).reset_index(drop=True)

#%%


#pick newest version of assembly
genomes=bdf["Genome assembly ID"]
ftpurl="https://ftp.ncbi.nlm.nih.gov/genomes/all/"

output_folder=os.getcwd()
genome_urls=[]

for i in genomes:
    acc=i.replace("_","").split(".")[0]
    genome_urls.append(ftpurl+"/".join([acc[i:i+3] for i in range(0, len(acc), 3)])+"/")

for url in genome_urls:

    #normally you would do this with ftp but the module is not really working at the moment
    t=requests.get(url).text
    versions=np.array([i.split(">")[-1] for i in t.split('/</a>')])#[::2]])
    newest_version=versions[np.argmax([int((i+".-1000").split(".",1)[-1].split("_")[0]) for i in versions])] 
    print(newest_version)
    
    full_url=url+newest_version+"/"+newest_version+"_protein.faa.gz"
    print(full_url)
    output_file=Path(output_folder,newest_version+"_protein.faa.gz")
    
    
    urlretrieve(full_url,output_file)
    #time.sleep(20)

extract(output_folder)