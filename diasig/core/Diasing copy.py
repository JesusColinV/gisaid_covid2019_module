import os
import re
from statistics import variance
import pandas as pd
import json
import logging
import datetime
from concurrent.futures import ThreadPoolExecutor
from ..utils.dictionaries import age_quinquennia

class Analyzer:
    def __init__(self, path:str, sep:str = "\t"):
        """_summary_
            La funcion de inicialización realiza todo el proceso del analisis
        Args:
            path (str): la dirección del archivo a analizar
            stp (str): el separador del archivo en caso de ser un formato no estandarizado

        Returns:
            df: un resumen del archivo obtenido despues del preprocesamiento
        """
        logging.basicConfig(filename="app.log", level="DEBUG")
        try:
            dictionaries = self.__load_dict__
            df = self.__load_file__(path,sep)
            df = self.__complete__(df)
            df = self.__normalize__(df,dictionaries)
            return (True,df,None)
        except Exception as ex:
            logging.error(ex)
            return (False,None,ex)
        
    def __load_file__(self, path:str, sep:str) -> pd.DataFrame:
        """_summary_
            Cargamos el archivo a analizar
        Args:
            path (str): la dirección del archivo a analizar
            stp (str): el separador del archivo en caso de ser un formato no estandarizado

        Returns:
            df (DataFrame): el archivo a orginal a analizar
        """
        
        reading = {".csv":pd.read_csv,".xlsx":pd.read_excel}
        
        for k,v in reading.items():
            if bool(re.search(k,path)):
                df = v(path)
                return df
        try:  
            df = pd.read_csv(path,sep=sep)
        except:
            assert "El formato del archivo o el separador no es reconocido, introduce sep = \"{tu separador}\" como parametro"
    
    def __load_dict__() -> dict:
        """_summary_
            Cargamos los catalogos con la logica del analisis
        Args:
            None

        Returns:
            dictionaries: catalogos con la logica del analisis
        """
        dictionaries = dict()
        #dictionaries[f'{age_quinquennia}'] = age_quinquennia
        dir = "../dictionaries"
        dictionaries = os.listdir(dir)
        with open(f"{dir}/{dictionaries}", 'r',encoding="utf-8") as f:
            dictionaries[f[:-5]] = json.load(f)
        return dictionaries
    
    def __complete__(self, df:pd.DataFrame) -> pd.DataFrame:
        """_summary_
            Cargamos los catalogos con la logica del analisis
        Args:
            df (DataFrame): el archivo a orginal a analizar

        Returns:
            dictionaries: catalogos con la logica del analisis
        """
        delay = {2020:0,2021:53,2022:106}
        df_cleaned = df[df['Collection date'].str.len()!=7]
        df_cleaned['Year'] = [datetime.datetime.strptime(i,'%Y-%m-%d').year for i in df_cleaned['Collection date']]
        df_cleaned['date'] = [datetime.datetime.strptime(i,'%Y-%m-%d').isocalendar() for i in df_cleaned['Collection date']]
        df_cleaned["week"] = [df_cleaned.date.iloc[i][1]  for i in range(len(df_cleaned))]
        df_cleaned["week_cont"] = [df_cleaned.date.iloc[i][1] + delay[df_cleaned.date.iloc[i][0]] for i in range(len(df_cleaned))]
        df_cleaned['state'] = [i.split('/')[2]for i in df_cleaned['Location']]
        #columnas que si se usan
        return df_cleaned
    
    def __normalize__(self,df:pd.DataFrame,dictionaries:dict) -> pd.DataFrame:
        dic_normalizer = {category : None for category in ("state","group_age","group_patient_status","age")}
        executor = ThreadPoolExecutor(max_workers = len(dic_normalizer))
        df['variant_type'] =                executor.submit(self.__linage_normalized__,df['Lineage'],dictionaries["variant_types"])
        df[['age','group_age']] =           executor.submit(self.__age_normalized__,df['Patient age'],dic_normalizer,dic_normalizer)
        df['group_status'] =                executor.submit(self.__state_group_patient_status__,df['Patient age'],dic_normalizer)
        df[['state_key','region_key']] =    executor.submit(self.__state_normalized__,df['state'],dic_normalizer)
        return df
    
    def __age_normalized__(self,df,dic,dic_q)-> list:
        """_summary_
            noramliza edad y edge group
        Args:
            df (_type_): _description_
            dic (_type_): _description_

        Returns:
            list: _description_
        """
        variant = list()
        for i in range(len(df)):
            variant.append(int(df.iloc[i])) 
            for k,v in dic.items():
                if re.search(df.iloc[i],v):
                    for q in dic_q.keys():
                        if k in dic_q[q]:
                            variant.append((k,q))
        return variant

    def __state_normalized__(self,df,dic) -> list:
        """_summary_
            clave y valor real del estado unificado
        Args:
            df (DataFrame): _description_
            dic (dict): _description_
        Returns:
            list<tuple>: clave y valor de estados 
        """     
        variant = list()
        for i in range(len(df)):
            k = variant[df['state'].iloc[i]]
            if k in dic.keys(): 
                variant.append((k,dic[k]))
            else:
                variant.append((99,dic[99]))
        return variant
        
    
    def __state_group_patient_status__(self,df:pd.DataFrame, dic:dict) -> list:
        variant = list()
        for i in range(len(df)):
            variant.append(dic[df['Patient status'].iloc[i].lower()])
        return variant
    
    def __linage_normalized__(self, df:pd.DataFrame, dic:dict) -> list:
        variant = list()
        for i in range(len(df)):
            if re.search('AY',df.iloc[i]):
                variant.append('Delta')
            else:
                for k,v in dic.items():
                    
                    if df.iloc[i] in v:
                        variant.append(k)
                        break
                    if k == 'Other linages':
                        variant.append(k)
                        break
        return variant