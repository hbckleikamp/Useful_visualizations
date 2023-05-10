# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:17:03 2023

@author: ZR48SA
"""

#Source: https://gist.github.com/tommycarstensen/ec3c57761f3846c339de925b66f4ac1b
#SOURCE: https://github.com/bibcure/title2bib/blob/master/title2bib/bin/title2bib
#Source: https://github.com/fphammerle/pubmed-bibtex/blob/9cd3a7c057dbddd74d6a2147451a70c2b7f0c350/pubmed_bibtex/__init__.py#L59


import requests
import numpy as np
import xml.etree.ElementTree as ET
import sys
import calendar
import io
import pandas as pd
from unidecode import unidecode




#%%

file="C:/MultiNovo/For publication/Selenium/References.txt" #file containing mla-style references
output_file="refs.bib"


bare_url = "http://api.crossref.org/"
lines=pd.read_csv(file,sep="\t",header=None)



titles=[]
ids=[]
bibs=[]
not_found=[]
c=0
for line in lines.loc[:,0]:
    c+=1
    print(c)
    sims=[]


    if '"' in line:
        title=line.split('"')[1].strip()
    else:
        tits=line.split(".") #longest string between two dots, strip space and quotes
        title=tits[np.argmax(np.array([len(i) for i in tits]))].strip(' "')
    

    #First try with cross-ref
    found=False
    params = {"query.bibliographic": line, "rows": 20}
    print(params)

    
    url = bare_url+"works"
    r= requests.get(url, params=params)
    items = r.json()["message"]["items"]
    
    
    
    for i, item in enumerate(items):
    
        # if c==4:
        #     "i"+1    
    
        if "title" in item.keys() and "author" in item.keys():
           
            
            title_item = unidecode(item["title"][0])
            try:
                
                title_item = title_item.decode("utf-8")
    
    
            except:
                pass
            
            
                
            item["title"] = title_item
            
            #word similarity based on title, authors, journal and year
            first_author=[str(i.get("family")).lower() for i in item["author"]] # if i.get("sequence")=="first"]
            year="("+str(item["issued"]["date-parts"][0][0])+")"
            test=set([i.strip('" .,:') for i in title_item.lower().split()]+first_author+[year]) 
            sline=set([i.strip('" .,:') for i in line.lower().split()])
            
            
            #wd=len(set(test)-set(sline)) #word distance
            #find entry with longest set inner
            wd = len(set(test) & set(sline))
   


            
            
            if title_item in title.lower(): #or wd==0:
                if year in line: #check if year is correct
                    bibs.append(item)
                    print(item["title"])
                    found=True
                    break
        else:
            wd=0
            
        sims.append(wd)
            
    if not found:

        bibs.append(items[np.argmax(np.array(sims))])
        #pick most similar entry

headers = {
    'Accept': 'application/x-bibtex;q=0.5',
}


#%%
btex=[]
for ix,bib in enumerate(bibs):
    print(ix)
    btex.append(requests.get(bib.get("URL"), headers=headers).text) 
    
#%%
with io.open(output_file, 'w', encoding = "utf-8") as bibfile:
    bibfile.write("\n".join(btex)+"\n")

         
        






  
            #%%    
    # if not found:

    #     #then try with entrez
        

    #     print(title)
    
    #     titles.append(title)
    
        
    #     url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&field=title&term="+title
    
    #     r=requests.get(url)
    #     pmids=[i.split("</Id>")[0] for i in r.text.split("<Id>")[1:]]

       

    #     # You can either get information with texmed api or entrez
    #     # for pmid in pmids:
    #     #     url="https://www.bioinformatics.org/texmed/cgi-bin/list.cgi?PMID="+str(pmid)+"&linkOut"
    #     #     r=requests.get(url)
        
        
        
    #     # #https://gist.github.com/tommycarstensen/ec3c57761f3846c339de925b66f4ac1b
    #     efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    #     r = requests.get(
    #         '{}?db=pubmed&id={}&rettype=abstract'.format(efetch, ','.join(pmids)))
            

    
    #     t=r.text
        
    #     ## Loop over the PubMed IDs and parse the XML.
    #     root = ET.fromstring(r.text)
    #     for PubmedArticle in root.iter('PubmedArticle'):
            
            
    #         ArticleTitle = PubmedArticle.find('./MedlineCitation/Article/ArticleTitle')
    #         if unidecode(ArticleTitle.text.lower()) in unidecode(line.lower()):
                
    #             PMID = PubmedArticle.find('./MedlineCitation/PMID')
    #             ISSN = PubmedArticle.find('./MedlineCitation/Article/Journal/ISSN')
    #             Volume = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/Volume')
    #             Issue = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/Issue')
    #             Year = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/PubDate/Year')
    #             Month = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/PubDate/Month')
    #         ##    Year = PubmedArticle.find('./MedlineCitation/Article/ArticleDate/Year')
    #         ##    Month = PubmedArticle.find('./MedlineCitation/Article/ArticleDate/Month')
        
    #             Title = PubmedArticle.find('./MedlineCitation/Article/Journal/Title')
    #             MedlinePgn = PubmedArticle.find('./MedlineCitation/Article/Pagination/MedlinePgn')
    #             Abstract = PubmedArticle.find('./MedlineCitation/Article/Abstract/AbstractText')
    #             authors = []
    #             for Author in PubmedArticle.iter('Author'):
    #                 try:
    #                     LastName = Author.find('LastName').text
    #                     ForeName = Author.find('ForeName').text
    #                 except AttributeError:  # e.g. CollectiveName
    #                     continue
    #                 authors.append('{}, {}'.format(LastName, ForeName))
    #             ## Use InvestigatorList instead of AuthorList
    #             if len(authors) == 0:
    #                 ## './MedlineCitation/Article/Journal/InvestigatorList'
    #                 for Investigator in PubmedArticle.iter('Investigator'):
    #                     try:
    #                         LastName = Investigator.find('LastName').text
    #                         ForeName = Investigator.find('ForeName').text
    #                     except AttributeError:  # e.g. CollectiveName
    #                         continue
    #                     authors.append('{}, {}'.format(LastName, ForeName))
    #             if Year is None:
    #                 _ = PubmedArticle.find('./MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate')
    #                 Year = _.text[:4]
    #                 Month = '{:02d}'.format(list(calendar.month_abbr).index(_.text[5:8]))
    #             else:
    #                 Year = Year.text
    #                 if Month is not None:
    #                     Month = Month.text
    #             try:
    #                 for _ in (PMID.text, Volume.text, Title.text, ArticleTitle.text, MedlinePgn.text, Abstract.text, ''.join(authors)):
    #             ##        assert '"' not in _, _
    #                     if _ is None:
    #                         continue
    #                     assert '{' not in _, _
    #                     assert '}' not in _, _
    #             except AttributeError:
    #                 pass
                
                
                
    #             ## Print the bibtex formatted output.
    #             try:
    #                 print('@Article{{{}{}pmid{},'.format(
    #                     authors[0].split(',')[0], Year, PMID.text))
    #             except IndexError:
    #                 print('IndexError', pmids, file=sys.stderr, flush=True)
    #             except AttributeError:
    #                 print('AttributeError', pmids, file=sys.stderr, flush=True)
    #             print(' author="{}",'.format(' AND '.join(authors)))
    #             print(' title={{{}}},'.format(ArticleTitle.text))
    #             print(' journal={{{}}},'.format(Title.text))
    #             print(' year={{{}}},'.format(Year))
    #             if Volume is not None:
    #                 print(' volume={{{}}},'.format(Volume.text))
    #             if Issue is not None:
    #                 print(' number={{{}}},'.format(Issue.text))
    #             if MedlinePgn is not None:
    #                 print(' pages={{{}}},'.format(MedlinePgn.text))
    #             if Month is not None:
    #                 print(' month={{{}}},'.format(Month))
    #         ##    print(' Abstract={{{}}},'.format(Abstract.text))
    #             print(' ISSN={{{}}},'.format(ISSN.text))
    #             print('}')
        
    #             found=True
                
    #             "i"+1
                


                
                

