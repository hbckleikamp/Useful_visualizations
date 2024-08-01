# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 11:27:46 2024

@author: hkleikamp
"""

import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))

print(os.getcwd())

#%%

import pandas as pd
import numpy as np
from io import StringIO
import os
from pathlib import Path

from urllib.request import urlretrieve
import gzip, zipfile, shutil, tarfile, requests, time



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
    
#read input file
file="C:/Users/hkleikamp/Downloads/SI3_List of organisms selected for building DB48 database.xlsx"
df=pd.read_excel(file,engine="openpyxl")
c=df.iloc[1,:]
df=df.iloc[2:]
df.columns=c
taxids=df["Species Taxid"].values.astype(str)

#Get proteomes


url_l="https://rest.uniprot.org/proteomes/stream?compressed=false&fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd%2Cgenome_assembly%2Cgenome_representation&format=tsv&query=%28%28taxonomy_id%3A"
url_r="%29%29"

uniparc_l="https://rest.uniprot.org/uniparc/stream?compressed=true&format=fasta&query=%28%28organism_id%3A"
uniparc_r="%29%29"



proteomes=[]
for ixt,taxid in enumerate(taxids[11:]):
    print(ixt)
    
    


    #%%    
    print(taxid)
    ps=[]
    for ix,i in enumerate(["reference","other","redundant","excluded"]):
        url=url_l+taxid+url_r+"AND+%28proteome_type%3A" +str(ix+1)+"%29"
      
        r=requests.get(url)
    
        proteomes=pd.read_csv(StringIO(r.text),sep="\t")
        
        if len(proteomes):
            proteomes["type_name"]=i
            proteomes["type"]=ix
    
        ps.append(proteomes)
        
    proteomes=pd.concat(ps)

    
    
    
    
    #%%
    
    output_folder=os.getcwd()
    
    
    if len(proteomes):
        
        #parse busco
        proteomes["BUSCO"]=proteomes["BUSCO"].str.replace("[",",",regex=False)
        proteomes["BUSCO"]=proteomes["BUSCO"].str.replace("]","",regex=False)
        busco=(proteomes["BUSCO"]+["%"]).str.rsplit(",",expand=True).applymap(lambda x: x[2:-1]).astype(float)
        busco.columns=["completeness","single","duplicated","fragmented","missing","genes"]
        
        #pick most complete proteome
        df=pd.concat([proteomes,busco],axis=1)
        
        if len(df[df["Organism Id"]==int(taxid)]):
            df=df[df["Organism Id"]==int(taxid)]
        
        df=df.sort_values(by=["type","completeness",'Protein count'],ascending=[True,False,False])
        bdf=df.iloc[0,:].T

    
    
        uniprot_url="https://rest.uniprot.org/uniprotkb/search?compressed=false&format=tsv&query=%28%28proteome%3A"+bdf["Proteome Id"]+"%29%29&size=2"
        
        if len(pd.read_csv(StringIO(requests.get(uniprot_url).text),sep="\t")):
            
            full_url="https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28proteome%3A"+bdf["Proteome Id"]+"%29%29"
            output_file=str(Path(output_folder,bdf["Proteome Id"]+"_UniProt_Proteome.faa.gz"))
            
            
        else:
            
            full_url="https://rest.uniprot.org/uniparc/stream?compressed=true&format=fasta&query=%28%28proteome%3A"+bdf["Proteome Id"]+"%29%29"
            output_file=str(Path(output_folder,bdf["Proteome Id"]+"_UniParc_Proteome.faa.gz"))
            
            
      
        
    

    else:
    
        full_url=uniparc_l+bdf["Proteome Id"]+uniparc_r
        output_file=str(Path(output_folder,taxid+"_UniParc_Proteins.faa.gz"))
    
    
    print(full_url)
    
    
    
    uncompressed=output_file.replace("faa.gz","faa")
    
    urlretrieve(full_url,output_file)
 
    extract(output_folder)

    
    #%% Parse database 
    

    Equate_IL=True         # change I and J into L 
    Remove_ambiguous=True  # remove ambiguous amino acids "B","X","Z","[","(" , and J in case IL is not equated
    No_Fragments=True     # remove incomplete sequences from UniprotKB contain (Fragment)/(Fragments) in header
    Add_decoy=False        # append decoy of reversed or scrambled peptides
    Add_taxid=True         # add taxonomy id to header description
    
    
    Ambiguous_AAs=["B","O","U","X","Z","[","("]
    decoy_delimiter="decoy_"
    decoy_method="reverse" #or "scramble"
    
    
    def chunk_gen(it,size=10**6):
        c=itertools.count()
        for _,g in itertools.groupby(it,lambda _:next(c)//size):
            yield g
        
    import Bio
    from Bio import SeqIO
    import itertools
    import random

    if any([Equate_IL,Remove_ambiguous,No_Fragments,Add_taxid]):

        Output_path=uncompressed
        suf=Path(Output_path).suffix
        if Remove_ambiguous: Output_path=Output_path.replace(suf,"_NoAmb"+suf)
        if No_Fragments:     Output_path=Output_path.replace(suf,"_NoFrag"+suf)
        if Equate_IL:        Output_path=Output_path.replace(suf,"_IJeqL"+suf)
        if Add_decoy:        Output_path=Output_path.replace(suf,"_Decoy"+suf)
        if Add_taxid:        Output_path=Output_path.replace(suf,"_taxid"+suf)        

        recs=SeqIO.parse(uncompressed,format="fasta")
        chunks=chunk_gen(recs)

   
        print("writing "+Path(Output_path).stem)
        with open(Output_path,"w+") as f:
            for c in chunks: 
       
                chunk_df=pd.DataFrame([[r.id,str(r.seq),r.description] for r in c],columns=["id","seq","description"])
                
                    
                
                
                if Add_taxid:        chunk_df.loc[~chunk_df["description"].str.contains("OX="),"description"]+=" OX="+str(taxid)    
                if Equate_IL:        chunk_df["seq"]=chunk_df["seq"].str.replace("I","L").str.replace("J","L")
                if Remove_ambiguous: chunk_df=chunk_df[~pd.concat([chunk_df["seq"].str.contains(aa,regex=False) for aa in Ambiguous_AAs],axis=1).any(axis=1)]
                if No_Fragments:     chunk_df=chunk_df[~chunk_df["description"].str.contains("(Fragment",regex=False)]
    
        
                if Add_decoy:
                    decoy=chunk_df.copy()
                    if decoy_method=="scramble": decoy["seq"]=decoy.seq.apply(lambda x: ''.join(random.sample(x,len(x))))
                    if decoy_method=="reverse":  decoy["seq"]=decoy.seq.apply(lambda x: x[::-1])
                    decoy["id"]=decoy_delimiter+decoy["id"]
                    chunk_df=pd.concat([chunk_df,decoy])
                
               
                f.write("\n"+"\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"]))
    
 