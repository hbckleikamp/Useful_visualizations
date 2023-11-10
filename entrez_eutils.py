# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 11:52:11 2023

@author: hkleikamp
"""


import requests


Term='"Methyl-coenzyme M reductase subunit alpha"'




#%% Get gene IDs from Term
retmax=1000
url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&retmax="+str(retmax)+"&term="+Term



r=requests.get(url).text
ids=[i.split("</Id>")[0] for i in r.split("<Id>")[1:]]
tc=r.split("<Count>")[1].split("</Count>")[0]


offset=0
for i in range(int(tc)//retmax):
    print(i)
    offset+=retmax
    nurl=url+"&retstart="+str(offset)
    r=requests.get(nurl).text
    ids+=[i.split("</Id>")[0] for i in r.split("<Id>")[1:]]

ids=list(set(ids))

#%% Get fastas from gene IDS (simple)



# efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
# batch=200
# chunks=[ids[i:i + batch] for i in range(0, len(ids), batch)]

# with open("mcra.fa","w") as f:
#     for ix,chunk in enumerate(chunks):
#         print(ix)        
#         url='{}?db=nuccore&id={}&rettype=fasta&retmode=text'.format(efetch, ','.join(chunk))
#         r = requests.get(url)
#         f.write(r.text)
        
#%% Get fastas and taxids and orgname

import pandas as pd

efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
batch=200
chunks=[ids[i:i + batch] for i in range(0, len(ids), batch)]

with open("mcra.fa","w") as f:

    for ix,chunk in enumerate(chunks):
        print(ix)
        url='{}?db=nuccore&id={}&rettype=fasta&retmode=xml'.format(efetch, ','.join(chunk))
        r = requests.get(url).text
        
        targets=["TSeq_accver","TSeq_defline","TSeq_orgname","TSeq_taxid","TSeq_sequence"]
        prefix=[">"," "," OS="," OX=","\n"]
        
        f.write("\n".join(pd.DataFrame([[prefix[tix]+i.split("</"+t+">")[0] for i in r.split("<"+t+">")[1:]] for tix,t in enumerate(targets)]).T.sum(axis=1).tolist()))
        
     
    
