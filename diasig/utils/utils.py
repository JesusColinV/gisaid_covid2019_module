import numpy as np
import re
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

class Counter:
    def __init__(self) -> None:
        pass
   
    @staticmethod
    def get_mutations(df,parts_protein:dict,protein:tuple = ("Spike_","E_","N_","M_","NSP")) -> dict:
        
        position = lambda val : int(re.findall("\d+",val.split('_')[1])[0])
        
        table = lambda val,count,aa: [(val[i],count[i],val[i].split('_')[0],position(val[i]), p_p(val,aa)) for i in range(len(val))]
        table2 = lambda val,count: [(val[i],count[i],val[i].split('_')[0],position(val[i])) for i in range(len(val))]
        v_c = lambda val: np.unique(val, return_counts=True)
        p_p = lambda val,aa:  [pp[0] for pp in parts_protein[aa] if position(val) in range(pp[1],pp[2])][0]
        
        with ThreadPoolExecutor(max_workers=len(protein)) as executor:
            futures = list()
            unique = list()
            for aa in protein[:-1]:
                futures.append(executor.submit(search_amino,aa,df['AA Substitutions']))
                values = [f.result() for f in futures]
            for i,aa in enumerate(protein[:-1]):
                v,c = v_c(values[i])
                unique.append(executor.submit(table,v,c,aa[:-1]))
            unique.append(executor.submit(table2,v,c))
        amino = {protein[i]:u.result() for i,u in enumerate(unique)}
        return amino
    
    @staticmethod
    def get_table(amino,parts_protein:dict,txt:str = ""):
        #, #Guruprasad L. Human SARS CoV-2 spike protein mutations. Proteins. 2021 May;89(5):569-576. doi: 10.1002/prot.26042. Epub 2021 Jan 17. PMID: 33423311; PMCID: PMC8014176
        futures = list()
        with ThreadPoolExecutor(max_workers=len(amino.keys())) as executor:
            cut = lambda tpl,i : [t[i] for t in tpl]
            for i in range(len(amino.keys())):
                f = [futures.append(executor.submit(cut,amino,j)) for j in range(5)]


        #    
        #segment,"count","full","change","position"
        #title = lambda segment,sep ="" : { t+"_"+sep:list() for t in (segment,"count","full","change","position")}
        #val = title("amigo","2")
        #myDF=[]
        #{}
        #for segment in amino.keys():
        #    myDF.append((amino[segment],[segment,'count'+txt,'full'+txt,'change'+txt,'position'+txt])) 
        #for i in range(len(myDF)-1):
        #    myDF[i]['area'] = [x[0] for num in myDF[i].position for x in parts_protein[myDF[i].full[0]] if num >= x[1] and num <= x[2]]
        #return myDF
    @staticmethod
    def properties_mutations(df,aminos:list):
        df_amino_changes = dict()
        for i in range(len(df)):
            #for segment in protein:
            for change in aminos:
                if change not in df_amino_changes.keys():
                    df_amino_changes[change] = {i:0 for i in df['variant_type'].unique()}       
                if re.search(change,df['AA Substitutions'].iloc[i]) is not None:
                    df_amino_changes[change][df['variant_type'].iloc[i]] +=1
        return df_amino_changes

def search_amino(segment:str,sustitution, amino:list = []) -> list:
        for user in sustitution:
            for aa in user.replace(r'(','').replace(r')','').split(','):
                if bool(re.search(segment,aa)):
                    amino.append(aa)
                    break
        return amino  