#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 14:51:52 2021

@author: hugokleikamp
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

os.chdir(Path(__file__).parents[0])
print(os.getcwd())

#%%
import pandas as pd    
import numpy as np  
  
#%% load to dataframes

ranks=["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
taxonomies=["gtdb_taxonomy","ncbi_taxonomy","ssu_silva_taxonomy"]


gtdbfiles=["/Volumes/Seagate_SSD/third_paper_newest/Data_resources/Metadata/GTDB/bac120_metadata_r202.tsv",
"/Volumes/Seagate_SSD/third_paper_newest/Data_resources/Metadata/GTDB/ar122_metadata_r202.tsv"]

mapdf = pd.DataFrame()
for file in  gtdbfiles:
    mapdf=mapdf.append(pd.read_csv(file, sep='\t', usecols=["accession"]+taxonomies))    
mapdf=mapdf[~(mapdf=="none").any(axis=1)]        
mapdf=mapdf.reset_index(drop=False)

#%% fix silva taxonomy to only contain only specified ranks with prefixes

stax_r=pd.read_csv(
"/Volumes/Seagate_SSD/third_paper_newest/Data_resources/Metadata/Silva/tax_slv_ssu_138.1.txt", sep='\t', header=None)
stax_r.columns=["silva_lineage","silva_taxid","silva_rank","taxon_lvl","version"]


def get_taxon_rank(x):
    sx=x.split(";")
    if len(sx[-1]):
        x=sx[-1]
    else:
        x=sx[-2]
    return x

stax_r["taxon_lvl"]=stax_r["silva_lineage"].apply(lambda x: get_taxon_rank(x))
stax_r["prefix"]=stax_r["silva_rank"].apply(lambda x: x[0]+"__")+stax_r["taxon_lvl"]

silva_tax=mapdf["ssu_silva_taxonomy"].str.rsplit(";",expand=True)

#%%
not_species=stax_r["taxon_lvl"].values.tolist()

species=silva_tax.mask(silva_tax.isin(not_species))


def get_species(x):
    #print(x)
    species_rank=[i for i in x if i==i and i != None]  #this assumes that the highest rank is species
    if len(species_rank):
        species_rank=species_rank[-1]
    else:
        species_rank=""

    return species_rank

species=species.apply(lambda x: get_species(x),axis=1)


silva_ranks=["domain", "phylum", "class", "order", "family", "genus", "species"]
stax_r=stax_r[stax_r["silva_rank"].isin(silva_ranks)]
allowed_taxa=stax_r["taxon_lvl"].values


silva_tax=silva_tax[silva_tax.isin(allowed_taxa)].reset_index(drop=True)
silva_tax=silva_tax.iloc[:,0:7]
silva_tax.columns=silva_ranks
silva_tax["species"]=species.tolist()


#index problems?
silva_tax=silva_tax.mask(silva_tax.isnull(),"")
silva_tax.replace(np.nan,"",inplace=True)
silva_tax=silva_tax.reset_index(drop=True)
#silva_tax=silva_tax.applymap(lambda x: "" if x.isnull())


rank_prefix=["d__", "p__", "c__", "o__", "f__", "g__", "s__"] #add ranks
for ix, i in enumerate(silva_tax.columns.tolist()):
    silva_tax[i]=rank_prefix[ix]+silva_tax[i]
    

mapdf["ssu_silva_taxonomy"]=silva_tax.apply(lambda x: ";".join(x),axis=1)
mapdf=mapdf.reset_index(drop=True)

#%% map silva and ncbi to GTDB      

#expand taxonomies

taxas=[]
for i in taxonomies: 
    taxa=[i+"_"+rank for rank in ranks]
    taxas.append(taxa)
    mapdf[taxa]=mapdf[i].str.rsplit(";",expand=True)

gtdbtax=taxas[0]
ncbitax=taxas[1]
silvatax=taxas[2]





#%% From GTDB 2 DB

def most_frequent_taxa_nonz(x):
    c=x.value_counts()
    if len(c):
        x=[c.index[0].strip(),c[0], c[0]/sum(c)]
        
        ind=0
        while ind<len(c) :
            if len(c.index[ind].strip())>3: #non p__ or "" including single whitespace 
                x=[c.index[ind].strip(), c[ind], c[ind]/sum(c)]
                break
            ind=ind+1
    else:
        x=["",0,0]
    return x

conv=pd.DataFrame()  
for ix,tax in enumerate([gtdbtax,ncbitax,silvatax]):
    alltax=pd.Series()
    alllin=pd.Series()
    prefix="_".join(tax[0].split("_")[0:-1])+"_"
    
    
    for iy,i in enumerate(gtdbtax):
        alltax=alltax.append(mapdf.groupby(i)[tax[iy]].apply(lambda x: most_frequent_taxa_nonz(x)))
    
        linranks=tax[0:iy+1]
        lin=mapdf[linranks].apply(lambda x: ";".join(x),axis=1)
        alllin=alllin.append(lin.groupby(mapdf[linranks[-1]]).apply(lambda x: most_frequent_taxa_nonz(x)))
        

    conv[[prefix+"tax",prefix+"taxcount",prefix+"taxfreq"]]=pd.DataFrame(alltax.values.tolist()
                                                                         ,columns=[prefix+"tax",prefix+"taxcount",prefix+"taxfreq"])
    conv=conv.merge(pd.DataFrame(alllin.values.tolist()
                                 ,columns=[prefix+"lin",prefix+"lincount",prefix+"linfreq"]).set_index(alllin.index).reset_index(drop=False).rename(columns={"index":prefix+"tax"}),on=prefix+"tax",how="left")

GTDB2db=conv
GTDB2db.to_excel("GTDB2dbs_mapping.xlsx")

#%% From db to GTDB




conv=pd.DataFrame()  
for ix,tax in enumerate([ncbitax,silvatax]):
    alltax=pd.Series()
    alllin=pd.Series()
    prefix="_".join(tax[0].split("_")[0:-1])+"_"

 
    
    for iy,i in enumerate(gtdbtax):
        alltax=alltax.append(mapdf.groupby(tax[iy])[gtdbtax[iy]].apply(lambda x: most_frequent_taxa_nonz(x)))
    
        linranks=tax[0:iy+1]
        lin=mapdf[linranks].apply(lambda x: ";".join(x),axis=1)
        alllin=alllin.append(lin.groupby(mapdf[linranks[-1]]).apply(lambda x: most_frequent_taxa_nonz(x)))

    
    at=pd.DataFrame(alltax.values.tolist())
    at.index=alltax.index
    
    at=at.reset_index()
    at.columns=[prefix+"_taxonomy_tax","gtdb_taxonomy_tax", "gtdb_taxonomy_taxcount", "gtdb_taxonomy_taxfreq"]
    at.to_excel(prefix+"_to_GTDB_mapping.xlsx")
