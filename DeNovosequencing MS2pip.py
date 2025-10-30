# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 11:29:37 2022

@author: ZR48SA
"""

#%% modules
import numpy as np
import pandas as pd
import itertools

#why not correct:
    #only up to for amino acid combinatnions (could have further N/GG Q/AG/GA substitutions)
    #doesnt consider modifications
    #ion pairs often are absent simultenaeously

#%% Define masses

mass_H=1.00782503223
mass_C=12
mass_N=14.00307400443
mass_O=15.99491461957
mass_P=30.97376199842
mass_S=31.9720711744
mass_Se=79.91652
ele_mass=[mass_H,mass_C,mass_O,mass_N,mass_P,mass_S,mass_Se]
elements=["Hydrogen","Carbon","Oxygen","Nitrogen","Phosphor","Sulfur","Selenium"]

aa_molecular_formulas=pd.DataFrame([
#AA  H    C O N P S Se 
#["", 0,   0,0,0,0,0,0],
["A",5,   3,1,1,0,0,0],
#["C",6,   3,1,1,0,1,0],
["C",7,5,2,2,0,1,0], #carbamidomethylated cysteine "H(3) C(2) N O"
["D",5,   4,3,1,0,0,0],
["E",7,   5,3,1,0,0,0],
["F",9,   9,1,1,0,0,0],
["G",3,   2,1,1,0,0,0],
["H",7,   6,1,3,0,0,0],
#["I",12,  6,1,1,0,0,0],
#["J",11,  6,1,1,0,0,0],
["K",12,  6,1,2,0,0,0],
["L",11,  6,1,1,0,0,0],
["M",9,   5,1,1,0,1,0],
["N",6,   4,2,2,0,0,0],
["P",7,   5,1,1,0,0,0],
["Q",8,   5,2,2,0,0,0],
["R",12,  6,1,4,0,0,0],
["S",5,   3,2,1,0,0,0],
["T",7,   4,2,1,0,0,0],
["V",9,   5,1,1,0,0,0],
["W",10, 11,1,2,0,0,0],
["Y",9,   9,2,1,0,0,0],
# ["U",5,   3,1,1,0,0,1],
# ["O",19, 12,2,3,0,0,0], #no side chain neutral losses are described for U and O in expert system
#modified residues
# ["Mox", 9,5,2,1,0,1,0], #oxidized methionine
# ["Ccam",8,5,2,2,0,1,0], #carbamidomethylated cysteine "H(3) C(2) N O"
# ["Sph", 6,3,5,1,1,0,0], #phosphorylation "H O(3) P"
# ["Tph", 8,4,5,1,1,0,0],
# ["Yph",10,9,5,1,1,0,0]],
],
    columns=["AA"]+elements).set_index("AA")
#add acetylations #no side chain neutral losses are described for acetylated amino acids
# ac=aa_molecular_formulas+[2,2,1,0,0,0,0] # "H(2) C(2) O"
# ac.index=aa_molecular_formulas.index+"ac"
# aa_molecular_formulas=pd.concat([aa_molecular_formulas,ac])
# aa_molecular_formulas=aa_molecular_formulas.T.to_dict(orient='list')

#fix carbamidomethyl
#aa_molecular_formulas.loc["C"]=aa_molecular_formulas.loc["Ccam"]
aa_mass=(aa_molecular_formulas*ele_mass).sum(axis=1)

#generate unique combinations with replacement (up to 4)
from itertools import combinations_with_replacement,permutations, product


combs=[]
for i in range(1,4): #not fully correct of course
    # combs.extend([c for c in combinations_with_replacement(range(len(aa_mass)),i)])
    # combs.extend([c for c in permutations(range(len(aa_mass)),i)])
    
    combs.extend(itertools.product(range(len(aa_mass)), repeat=i))
    
combs=pd.Series(combs).to_frame(name="aa").reset_index()
ce=combs.explode("aa")
ce["an"]=aa_mass.index[ce.aa.values.astype(int)]
ce["mass"]=aa_mass.values[ce.aa.values.astype(int)]

cem=ce.groupby("index")["mass"].sum().to_frame(name="mass")
cem["aa"]=ce.groupby("index")["an"].apply(list)
cem=cem.sort_values(by="mass").reset_index(drop=True)
um,uc=np.unique(cem.mass,return_counts=True)
cem.index=np.repeat(np.arange(len(uc)),uc)

cem.aa=cem.aa.apply("".join)

#unique permutations of a string instead of combinations with replacement.
#add number of unique combinations to column

from sklearn.neighbors import KDTree

tree = KDTree(um.reshape(1,-1).T, leaf_size=200) 

#vectorized find nearest mass
#https://stackoverflow.com/questions/8914491/finding-the-nearest-value-and-return-the-index-of-array-in-python
def find_closest(A, target): #returns index of closest array of A within target
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx


def zero_runs(a):
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

#chatgpt code
from itertools import combinations
  
def non_overlapping_full_coverage(ranges, gap_start, gap_end):
    """
    Find all minimal, non-overlapping combinations of ranges
    that fully cover the gap [gap_start, gap_end].
    """
    gap = set(range(gap_start, gap_end + 1))
    results = []
  
    def overlaps(r1, r2):
        a, b = r1
        c, d = r2
        return not (b < c or d < a)  # True if they overlap
  
    for r in range(1, len(ranges) + 1):
        for combo in combinations(ranges, r):
            # skip if any two ranges overlap
            if any(overlaps(combo[i], combo[j]) for i in range(len(combo)) for j in range(i+1, len(combo))):
                continue
  
            # compute total coverage
            covered = set()
            for a, b in combo:
                covered |= set(range(a, b + 1))
            if not gap.issubset(covered):
                continue
  
            # check minimality
            is_minimal = True
            for i in range(len(combo)):
                subcombo = combo[:i] + combo[i+1:]
                sub_covered = set()
                for x, y in subcombo:
                    sub_covered |= set(range(x, y + 1))
                if gap.issubset(sub_covered):
                    is_minimal = False
                    break
  
            if is_minimal:
                results.append(combo)
    return results

                    
#%%

filenames=["ecoli","yeast","pasteurella"]

msps=["C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/proteomes/ecoli_ms2pip_predictions.csv",
"C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/proteomes/yeast_ms2pip_predictions.csv",
"C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/proteomes/pasteurella_ms2pip_predictions.csv"]

for ixf,msp in enumerate(msps):
    
    cd=pd.read_csv(msp)
    cd["l"]=cd.spec_id.apply(len)-2
    cd["rion"]=cd["ionnumber"]
    cd.loc[cd.ion=="Y","rion"]=cd.loc[cd.ion=="Y","l"]-cd.loc[cd.ion=="Y","rion"]
    
    
    
    
    
    
    # Completeness
    hmass=1.007276
    water_mass=mass_H*2+mass_O
    ppm=20
    
    completeness=[1,0.9,0.8,0.7,0.6,0.5]
    
    counter=0
    
    expanded_peps=[]
    #%%
    sampleset=100 #2500 #larger space? 1000,10000
    if sampleset:
        cd=cd[cd.spec_id.isin(cd.spec_id.drop_duplicates().sample(sampleset))]
  
    for n,g in cd.groupby("spec_id"):
        

        #if n!='ADVQEK_2': continue
        #if n!="AAQMAAAGQTAQND_2": continue
        #if n!='AAAADTLLILGDSLSAGYR_2': continue #test
        # if n!="YWGQYSLFAQVLVGLLLANLALVLSR": continue 
        # else: break

        g=g.drop_duplicates()
        pre=g.groupby("rion")["mz"].sum()-hmass 
        
    
        gl=len(g)
        pep=g.spec_id.values[0].split("_")[0].replace("I","L")
        
        #if pep!='ACHWTVVLCFFLVALSGLSFFFPTLQWLTQTFGTPQMGR': continue
        
        #ms2pip mzs are actually incorrect??? (calculate correct masses)
        ladder=aa_mass.loc[list(pep)]
        ladder_B=np.cumsum(ladder.values)+hmass
        ladder_B[-1]+=water_mass
        ladder_Y=ladder.values[::-1]
        ladder_Y[0]+=(water_mass+hmass)
        ladder_Y=np.cumsum(ladder_Y)
        g=g.sort_values(by=["ion","ionnumber"])
        g["mz"]=np.hstack([ladder_B[:-1],ladder_Y[:-1]])
        lb=np.hstack([0,ladder_B])
        # lb=pd.DataFrame(np.vstack([np.arange(len(pep)).astype(int),ladder_B]).T,columns=["pos","m"]).set_index("m")
        # ly=pd.DataFrame(np.vstack([np.arange(len(pep)).astype(int),ladder_Y]).T,columns=["pos","m"]).set_index("m")
        
        g=g.sort_values(by="prediction",ascending=False)
        for x in completeness:
       
      
            g=g.iloc[:int(gl*x)].sort_values(by="ionnumber")
        
            #B ladder (partial)
            gb=g[g.ion=="B"]
            gbi=np.hstack([gb.ionnumber.values,len(pep)])-1
            cladder_B=np.diff(np.hstack([0,gb.mz,pre.values[0]]))
            cladder_B[0]-=hmass
            cladder_B[-1]-=water_mass #subtract water mass
        
            t=tree.query_radius(cladder_B.reshape(1,-1).T,r=cladder_B*ppm/1e6)
            tix=np.vstack([np.hstack([np.repeat(ix,len(i)) for ix,i in enumerate(t)]),np.hstack(t)]).T
            tdf=pd.DataFrame(tix,columns=["mx","cx"])
            tdf["pos"]=gbi[tdf.mx]
            
            solb=pd.DataFrame([np.repeat(tdf.pos.values,uc[tdf.cx]),cem.loc[tdf["cx"]].aa.values]).T
            solb.columns=["pos","aa"]
            solb["m"]=ladder_B[solb.pos.astype(int).values]
            solb["ml"]=solb["m"]-cem.loc[tdf.cx].mass.values #-hmass #mass span
            solb["posl"]=find_closest(lb,solb["ml"]) #position span (incorrect?)
            solb["c"]=solb.groupby("posl").transform("size")
            solb["i"]="b"
            
            #Y ladder (partial)
            gy=g[g.ion=="Y"]
            gyi=np.hstack([len(pep),gy.rion.values,])[::-1]-1
            cladder_Y=np.diff(np.hstack([0,g[g.ion=="Y"].mz,pre.values[0]]))[::-1]
            cladder_Y[-1]-=(water_mass+hmass) #subtract water mass & H+
          
            t=tree.query_radius(cladder_Y.reshape(1,-1).T,r=cladder_Y*ppm/1e6)
            tix=np.vstack([np.hstack([np.repeat(ix,len(i)) for ix,i in enumerate(t)]),np.hstack(t)]).T
            tdf=pd.DataFrame(tix,columns=["mx","cx"])
            tdf["pos"]=gyi[tdf.mx]
        
            soly=pd.DataFrame([np.repeat(tdf.pos.values,uc[tdf.cx]),cem.loc[tdf["cx"]].aa.values]).T
            soly.columns=["pos","aa"]
            soly["m"]=ladder_B[soly.pos.astype(int).values]
            
            soly["ml"]=soly["m"]-cem.loc[tdf.cx].mass.values #-hmass #mass span
            #soly["massl"]=soly["m"]-cem.loc[tdf.cx].mass.values
            soly["posl"]=find_closest(lb,soly["ml"]) #position span
            soly["c"]=soly.groupby("posl").transform("size")
            soly["i"]="y"
            
            solm=pd.concat([soly,solb])

            
            #check if every position can be solved:
            solm["span"]=(solm["pos"]-solm["posl"]+1).apply(np.arange)
            e=solm[["posl","span"]].explode("span")
            pepcombs=["unresolved"]
            no_combs=0 #placeholder
            
            #missing ions!
            if len(set(np.arange(len(pep)))-set((e.posl+e.span))): 
         
                #try final b/y diff solve
                l,r=solb.m.values[-1],soly.ml.values[0]
                if r>l:
                    d=np.array(r-l)
                    t=tree.query_radius(np.array(r-l).reshape(1,-1).T,r=d*ppm/1e6)
                    solby=cem.loc[np.hstack(t)]
                    base=solb.iloc[[-1]]
                    base.loc[:,["m","ml"]]+=solby.mass.values[0]
                    base.loc[:,["pos","posl"]]+=1
                    base.loc[:,"aa"]=solby.aa.values[0]
                    solm=pd.concat([solm,base])
                    solm["span"]=(solm["pos"]-solm["posl"]+1).apply(np.arange)
                    e=solm[["posl","span"]].explode("span")
                    
            solm["posr"]=len(pep)-solm["pos"]
            
            
            #check span completeness
            if not len(set(np.arange(len(pep)))-set((e.posl+e.span))): 
    
             
    
                usol=solm[["pos","aa","posl","c"]].drop_duplicates().reset_index(drop=True)
                zm=["*"]*len(pep) #left and right tag solve
                
                #fill_anchors
                for n,gr in usol.groupby("pos"):
                    sing=gr[gr.c==1]
                    #sing is not single!
                    
                    if len(sing)==1: 
                        for ia,a in enumerate(sing.aa.values):
                            zm[n+ia]=a
            
                #filter solutions recursive with borders (could also use regex or alignment?)
                exclude=[""]   #initialize
                while len(exclude): 
            
                    exclude=[]      
                    for n,gr in usol.groupby(["posl","pos"]):
                        rg="".join(zm[n[0]:n[1]+1])
                        if rg=="*": continue
                        aas=gr.aa
                        
                        ll,lr=0,0
                        if len(rg)==1:
                            b=aas.str.startswith(rg)
                        else:
                          
                            if "*" in rg:
                                srg=rg.split("*")
                                l,r=srg[0],srg[-1]
                                ll,lr=len(l),len(r)
                                b=aas.str.startswith(l) & aas.str.endswith(r)
                                
                                keep=gr[b]
                                if len(keep)==1:
                                    #get location of *
                                    if lr: zm[keep.posl.values[0]+ll]=keep.aa.values[0][ll:-(lr)]
                                    else:  zm[keep.posl.values[0]+ll]=keep.aa.values[0][ll:]
                            else:
                                b=aas.str.startswith(rg)
                    
                        b[~b].index.tolist()
                        exclude.append(b[~b].index.tolist())
                    exclude=np.hstack(exclude)    

                    #print(len(exclude))    
                    usol=usol[~usol.index.isin(exclude)]
                    
                #fill zm
            
                s=usol.groupby("pos").size()
                s1=usol[usol.pos.isin(s[s==1].index)]
                s1=s1[s1["pos"]==s1["posl"]]
                 
                for s in s1[["pos","aa"]].values:
                    zm[s[0]]=s[1]
                    
                ""+1
          
                # final zm cleanup
                if "*" in zm:
                    
                    a=np.hstack([(np.array(zm)!="*").astype(int),1]) #pad 1
                    
                    gap_info=zero_runs(a)
                    gd=np.diff(gap_info,axis=1).flatten()
                    
                    nonsec_gaps=gap_info[gd==1,0]
                    consec_gaps=gap_info[gd>1]
                    
           
            
                    for ngap in nonsec_gaps:
                        zm[ngap]=usol[usol.pos==ngap].aa.tolist()
                         

                    if len(consec_gaps):
                    
                        cgap_sols=[]
                        
                        for cgap in consec_gaps:
                            
                            gap_start, gap_end=cgap[0],cgap[1]-1
                            covers=usol[(usol.pos>=gap_start) & (usol.posl<=gap_end)]
                            cs=covers.groupby(["posl","pos"]).size()
                            ranges=covers[["posl","pos"]].drop_duplicates().values.astype(int).tolist()
                            pairs=non_overlapping_full_coverage(cs.index, gap_start, gap_end)
                            
                            
                            css=[np.prod(cs.loc[list(p)])  for p in pairs]
                            minimal_pairs=[i for ix,i in enumerate(pairs) if css[ix]==min(css)]
                            
                            #trim minimal pairs if they overlap
                            
                            sols=[]
                            
                            for mp in minimal_pairs:
                                
                                for pt in mp:
                                    
                                    p=covers.set_index(["posl","pos"]).loc[pt].reset_index()
                                    
                                    #left_trim,right_trim=zm[pt[0]:gap_start],zm[gap_end+1:pt[1]+1]
                                    left_trim,right_trim=gap_start-pt[0],pt[1]-gap_end
                                    
                                    if left_trim>0:
                                    
                                        if right_trim: p.aa=p.aa.apply(lambda x: x[left_trim:-right_trim])
                                        else:          p.aa=p.aa.apply(lambda x: x[left_trim:])
                                    sols.append(p.aa.tolist())
                                    
                                    #trimming not fully correct?
                               
                            cgap_sols.append(sols)    
                                
                        
     
                            fzm=[]
                            c=0
                            for ix,i in enumerate(zm):
                               if ix in consec_gaps[:,0]:
                                   fzm.append(["".join(i) for i in itertools.product(*sols)])
                                   c+=1
                               else:
                                   if i!="*":
                                       fzm.append(zm[ix])
                                       
                        else:
                            fzm=zm
                                   
                    no_combs=np.prod([len(i) for i in fzm])
                    if no_combs<100:
                        pepcombs=["".join(i) for i in itertools.product(*fzm)]
                                   
                

                else:        
                    pepcombs="".join(zm)




            
                #very weird buuuuggg
                if pep not in pepcombs:
                    ""+1 #test



#%%
                
#                     pepcombs=["".join(i) for i in itertools.product(*vc.groupby("posl")["aa"].apply(list).tolist())]
#                 else:
#                     pepcombs=["unspecific"]
                
#             if pep not in pepcombs:
#                 if pepcombs!=["unspecific"]:
#                     if pepcombs!=["unresolved"]:
#                         ""+1
            
#             #if len(exclude): ""+1
            
#                 #startswith, endswith, regex?
#                 #or alignment?
                
#                 # ""+1
#                 # g
            
#                 #fill zm with 
            
#             #%%
            
            
#             # exlude=[]
#             # for n,gr in sole.groupby("pos"):
#             #     sing=gr[gr.c==1]
                
 
#             #     #check previous span untill 0
                
                
                
#             #     if len(sing): 
                    
#             #         for ia,a in enumerate(sing.aa.values):
#             #             zm[n+ia]=a
                    
                
                
#                     #%%
                    
#                 #     lr.append(gr[gr.aa.str.startswith(sing.aa.values[0])])
                
                
                
                
#                 # else:         lr.append(gr)
#             lr=pd.concat(lr)
            
#             lr["c"]=lr.groupby(["posl","i"]).transform("size")
#             rl=[]
#             for n,gr in lr.groupby("posr"):
#                 sing=gr[gr.c==1]
#                 if len(sing): rl.append(gr[gr.aa.str.endswith(sing.aa.values[0])])
#                 else:         rl.append(gr)
#             rl=pd.concat(rl)
#             rl=rl[["pos","posl","aa","c"]].drop_duplicates()
       
#             d=(rl.pos-rl.posl+1)
#             rl["s"]=rl["c"]/d.values
#             rl["r"]=d.apply(np.arange)
#             rl=rl.explode("r")
#             rl["pos"]-=rl["r"]

#             #rl=rl.sort_values(by=["pos","s"]).groupby("pos",sort=False).nth(0)
#             rl=rl[(rl["s"]/rl.groupby("pos")["s"].transform("min"))==1] 
     
#             #vc=rl.groupby(rl.index).head(1)[["pos","aa"]].sort_values(by="pos") #not working? some index error
#             #vc=rl.groupby("pos").head(1)[["pos","aa"]].sort_values(by="pos") #not working? some index error  
#             vc=rl[["posl","aa"]].drop_duplicates().sort_values(by="posl") 
            
            
#             no_combs=np.prod(vc.groupby("posl").size())
           
            
#             if no_combs<100:
#                 pepcombs=["".join(i) for i in itertools.product(*vc.groupby("posl")["aa"].apply(list).tolist())]
#             else:
#                 pepcombs=["unspecific"]
                
#             if pep not in pepcombs:
#                 if pepcombs!="unspecific":
#                     ""+1
#             #%%
            
#             expanded_peps.append(np.hstack([np.tile(np.array([[counter,pep,x,no_combs]]),[len(pepcombs),1]),np.array(pepcombs).reshape(-1,1)]))
            
#         counter+=1
#         print(counter)
    
#     expanded_peps=pd.DataFrame(np.vstack(expanded_peps),columns=["ix","original","completeness","combinations","peptide"])
    
#     # expanded_peps.to_csv(filenames[ixf]+"_pepcombs.csv")
#     #submit expanded peps to uniPept
    
# #%%
#     peps=expanded_peps.peptide #[:100] #test
#     peps=peps[peps!="unspecific"]
    
#     upeps=expanded_peps.loc[expanded_peps.peptide=="unspecific","original"].drop_duplicates().to_frame()
#     upeps["taxon_rank"]="no rank"
#     upeps.columns=["peptide","taxon_rank"]
#     upeps["original"]=upeps["peptide"]
    
#     import random, re, requests
#     import threading, time, string
    
#     # Threading unipept
#     def unipept_scrape(r,url):
        
#        # r.append(pd.DataFrame(requests.get(url,stream=True).json())[["peptide","taxon_rank"]])
        
        
#         while True:
#             try:
#                 #r.extend(requests.get(url,stream=True).json())
                
#                 d=pd.DataFrame(requests.get(url,stream=True).json())
#                 if len(d): r.append(d[["peptide","taxon_rank"]])
#                 break
#             except:
#                 print("sleeping")
#                 time.sleep(2)
                
                
#     # chunker
#     def chunks(lst,n):
#         for i in range(0,len(lst),n):
#             yield lst[i:i+n]
    
    
#     twl='http://api.unipept.ugent.be/api/v1/pept2lca.json?input[]=';
#     #twr='&equate_il=true&extra=true&names=true';
#     twr='&equate_il=true&extra=true' #'&names=true';
    
#     comp_ranks=["domain","phylum","class","order","family","genus","species"]
#     ranks=[rank+"_name" for rank in comp_ranks]
#     fields=["peptide"]+["taxon_name"]+ranks
#     batchsize=100
    
#     #submit both to Unipept
    
    
    
#     for rev in range(2):
        
#         #reverse sequences and submit both to Unipept
#         if rev: peps=peps.str[:-1].str[::-1]+peps.str[-1] #reverse infront of C-term
#         #nontryp=((peps.str[-1]!="K") & (peps.str[-1]!="R")).sum()
#         base_thread=threading.active_count()
#         unipeps=peps
     
#         steps=list(range(0,len(unipeps),batchsize))
#         taxdf=[]
        
        
#         threads=[]
#         counter=0
#         for chunk in chunks(unipeps,batchsize):
#             ""+1
#             counter+=1
#             print(counter)
#             query="&input[]=".join(chunk)
#             turl=twl+query+twr 
   
#             time.sleep(0.1)
#             #without threading
#             #taxdf.append(pd.DataFrame(requests.get(turl,stream=True).json())[["peptide","taxon_rank"]])
            
#             #taxonomy
#             turl=twl+query+twr 
#             t=threading.Thread(target=unipept_scrape, args=[taxdf,turl])
#             t.start()
#             threads.append(t)
            
      
            
#             if not counter%20:
#                 print("unwinding, query at: "+str(round((counter*batchsize)/len(unipeps)*100,2))+"%")
#                 for thread in threads:
#                     thread.join()
#                 threads=[] #this seems to act different on windows?
             
#         for thread in threads:
#             thread.join()
                        
#         taxdf=pd.concat(taxdf)
#         taxdf=pd.concat([taxdf,upeps]).reset_index(drop=True)

        
        
#         # if rev: taxdf.to_csv(filenames[ixf]+"incomplete_reverse_unipept.csv")
#         # else:   taxdf.to_csv(filenames[ixf]+"incomplete_forward_unipept.csv")
    


# #%%

# names=["ecoli","yeast","pasteurella"]

# combs=["C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/ecoli_pepcombs.csv",
# "C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/yeast_pepcombs.csv",
# "C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/pasteurella_pepcombs.csv"]

# fow=[
# "C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/ecoliincomplete_forward_unipept.csv",
# "C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/yeastincomplete_forward_unipept.csv",
# "C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/pasteurellaincomplete_forward_unipept.csv"]


# rev=["C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/ecoliincomplete_reverse_unipept.csv",
# "C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/yeastincomplete_reverse_unipept.csv",
# "C:/Users/e_kle/Desktop/Antwerpen/FWO 2026/Case study/pasteurellaincomplete_reverse_unipept.csv"]


# rank_order=pd.Series([
# "no rank",
# "domain",
# "realm",
# "kingdom",
# "subkingdom",
# "superphylum",
# "phylum",
# "subphylum",
# "superclass",
# "class",
# "subclass",
# "infraclass",
# "superorder",
# "order",
# "suborder",
# "infraorder",
# "parvorder",
# "superfamily",
# "family",
# "subfamily",
# "tribe",
# "subtribe",
# "genus",
# "subgenus",
# "species_group",
# "species group",
# "species_subgroup",
# 'species subgroup',
# "species",
# "subspecies",
# "strain",
# "varietas",
# "forma"])

# #forward
# for ix,i in enumerate(combs):
    
#     cdf=pd.read_csv(combs[ix])
#     f,r=pd.read_csv(fow[ix]),pd.read_csv(rev[ix])
    
#     fm=cdf.merge(f[["peptide","taxon_rank"]],how="left",on="peptide").dropna()
#     fm["tax_c"]=rank_order.to_frame("r").reset_index().set_index("r").loc[fm.taxon_rank].values
    
#     fm=fm.sort_values(by=["ix","completeness","tax_c"])
#     fm=fm.groupby(["ix","completeness"],sort=False).nth(0)
#     fm["l"]=fm.original.apply(len)
    
    
#     #completeness
#     import matplotlib.pyplot as plt
#     fig,ax=plt.subplots()
#     for n,g in fm.groupby("completeness"):
#         s=pd.DataFrame(g.groupby("taxon_rank").size())
#         missin_taxa=list(set(comp_ranks)-set(s.index))
#         if missin_taxa:
#             s.loc[*list(set(comp_ranks)-set(s.index)),:]=0
        
#         rs=rank_order[rank_order.isin(s.index)].tolist()
#         s=s.iloc[:,0]
#         s=s[rs[::-1]].cumsum()[::-1]
#         d=s.loc[comp_ranks]
#         d=d/len(g)
#         plt.plot(np.arange(len(d)),d.values,label=n)

#     plt.legend()
#     ax.set_xticks(np.arange(len(d))) 
#     ax.set_xticklabels(d.index)
#     plt.title(names[ix])
#     plt.ylabel("Annotated fraction")
#     plt.savefig(names[ix]+"forward_completeness.png",dpi=300)
    
    


# #%% False positives

# fig,ax=plt.subplots()
# for ix,i in enumerate(combs):
    
#     cdf=pd.read_csv(combs[ix])
#     f,r=pd.read_csv(fow[ix]),pd.read_csv(rev[ix])
    
#     fm=fm.sort_values(by=["ix","completeness","tax_c"])
#     sf=fm[["ix","completeness"]].drop_duplicates().groupby("completeness").size()
    
#     rm=cdf.merge(r[["peptide","taxon_rank"]],how="left",on="peptide").dropna()
#     sr=rm[["ix","completeness"]].drop_duplicates().groupby("completeness").size()
    
#     s=sr/(sr+sf)

#     plt.scatter(s.index,s.values,label=names[ix]) 
#     plt.plot(s.index,s.values)  
    
# plt.ylabel("False positive hits")
# plt.xlabel("fragment completeness")
# plt.legend()
# plt.savefig("completeness_FDR.png",dpi=300)

# #%% FDR

# fig,ax=plt.subplots()
# for ix,i in enumerate(combs):
    
#     cdf=pd.read_csv(combs[ix])
#     f,r=pd.read_csv(fow[ix]),pd.read_csv(rev[ix])
    
#     rm=cdf.merge(r[["peptide","taxon_rank"]],how="left",on="peptide").dropna()
#     s=rm[["ix","completeness"]].drop_duplicates().groupby("completeness").size()

#     plt.scatter(s.index,s.values,label=names[ix]) 
#     plt.plot(s.index,s.values)  
    
# plt.ylabel("False positive hits")
# plt.xlabel("fragment completeness")
# plt.legend()
# plt.savefig("completeness_FP.png",dpi=300)


# #%% combs

# for ix,i in enumerate(combs):
    
#     cdf=pd.read_csv(combs[ix])
    
#     u=cdf[["original","completeness","combinations"]].drop_duplicates()
#     u["l"]=u["original"].apply(len)
    
#     fig,ax=plt.subplots()
#     for n,g in u.groupby("completeness"):
        
        
#         m=g.groupby("l")["combinations"].mean()
#         m=m[m.index<30]
#         plt.plot(m.index,m.values,label=n)
#         # plt.xlim(0,30)
#         # plt.ylim(0,100)
        
#     plt.legend()
    
    
# #%% probably hist is better

# for n,g in u.groupby("completeness"):
#     fig,ax=plt.subplots()
#     plt.hist(g.combinations,range=[0,20],bins=10)
#     plt.title(n)    
# #%% number of combinations
    
    