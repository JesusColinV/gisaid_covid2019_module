import os
import re
import pandas as pd
import json
import logging
import datetime
from concurrent.futures import ThreadPoolExecutor
from ..utils.dictionaries import age_quinquennia
from ..utils.utils import Counter

class Analyzer:
    def __init__(self, path:str, sep:str = "\t") -> None:
        """_summary_
            La funcion de inicialización realiza todo el proceso del analisis
        Args:
            path (str): la dirección del archivo a analizar
            stp (str): el separador del archivo en caso de ser un formato no estandarizado

        Returns:
            df: un resumen del archivo obtenido despues del preprocesamiento
        """
        self.path = path
        self.__name_df_cleaned = "filtered.csv"
        self.__name_df_completed = "completed.csv"
            
    def load(self) -> pd.DataFrame:
        logging.basicConfig(filename="app.log", level="DEBUG")
        try:
            dictionaries = self.__load_dict__()
            if self.__file_validator__(self.__name_df_completed):
                if self.__file_validator__(self.__name_df_cleaned):
                    df = self.__read_cleaned_file__()
                else:
                    df = self.__load_file__()
                    df = self.__complete__(df)
                    df.to_csv(self.__name_df_cleaned)
                df = self.__normalize__(df,dictionaries)
                df.to_csv(self.__name_df_completed)
                self.__counter__(df,dictionaries["proteins"])
        except Exception as ex:
            print(ex)
            logging.error(ex)
        return df
        
    def __counter__(self,df,dic):         
        amino_unique = Counter.get_mutations(df,dic)
        mutations = Counter.get_table(amino_unique,dic)
    def __load_file__(self) -> pd.DataFrame:
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
            if bool(re.search(k,self.path)):
                df = v(self.path)
                return df
        try:  
            df = pd.read_csv(self.path,sep="\t")
            return df
        except:
            assert "El formato del archivo o el separador no es reconocido, introduce sep = \"{tu separador}\" como parametro"
    
    def __file_validator__(self,file):
        files = os.listdir()
        return file in files
    def __read_cleaned_file__(self):
        return pd.read_csv(self.__name_df_cleaned)
    
    def __load_dict__(self) -> dict:
        """_summary_
            Cargamos los catalogos con la logica del analisis
        Args:
            None

        Returns:
            dictionaries: catalogos con la logica del analisis
        """
        dictionaries = dict()
        #dictionaries[f'{age_quinquennia}'] = age_quinquennia
        dir = "./diasig/dictionaries"
        files = os.listdir(dir)
        for file in files:
            with open(f"{dir}/{file}", 'r',encoding="utf-8") as f:
                dictionaries[file[:-5]] = json.load(f)
        return dictionaries
    
    def __complete__(self, df:pd.DataFrame) -> pd.DataFrame:
        """_summary_
            Cargamos los catalogos con la logica del analisis
        Args:
            df (DataFrame): el archivo a orginal a analizar

        Returns:
            dictionaries: catalogos con la logica del analisis
        """
        df_cleaned = df[df['Collection date'].str.len()>7]
        df_cleaned['date'] = [datetime.datetime.strptime(i,'%Y-%m-%d').isocalendar() for i in df_cleaned['Collection date']]
        with ThreadPoolExecutor(max_workers = 4) as executor:
            f1 = executor.submit(self.__year__,df_cleaned)
            f2 = executor.submit(self.__week__,df_cleaned)
            f3 = executor.submit(self.__week_cont__,df_cleaned)
            f4 = executor.submit(self.__state__,df_cleaned)
        df_cleaned['Year'] = f1.result()
        df_cleaned["week"] = f2.result() 
        df_cleaned["week_cont"] = f3.result() 
        df_cleaned['state'] = f4.result()
        return df_cleaned
    
    def __year__(self, df_cleaned:pd.DataFrame) -> list:
        return [datetime.datetime.strptime(i,'%Y-%m-%d').year for i in df_cleaned['Collection date']]
    
    def __week__(self, df_cleaned:pd.DataFrame) -> list:
        return [df_cleaned.date.iloc[i][1]  for i in range(len(df_cleaned))]
    
    def __week_cont__(self, df_cleaned:pd.DataFrame) -> list:
        delay = {2020:0,2021:53,2022:106}
        return [df_cleaned.date.iloc[i][1] + delay[df_cleaned.date.iloc[i][0]] for i in range(len(df_cleaned))]
    
    def __state__(self, df_cleaned:pd.DataFrame) -> list:
        return [i.split('/')[2]for i in df_cleaned['Location']]
    
    def __normalize__(self,df:pd.DataFrame,dictionaries:dict) -> pd.DataFrame:
        #dic_normalizer = {category : None for category in ("state","group_age","group_patient_status","age")}
        with ThreadPoolExecutor(max_workers = 4) as executor:
            f1 = executor.submit(self.__linage_normalized__,df['Lineage'],dictionaries["variant_types"])
            f2 = executor.submit(self.__age_normalized__,df[['Patient age','Gender']],dictionaries["age_unification"])
            f3 = executor.submit(self.__state_group_patient_status__,df['Patient status'],dictionaries["patient_status"])
            f4 = executor.submit(self.__state_normalized__,df['state'],dictionaries["states_types"],dictionaries["unique_states_types"])

        df['variant_type'] = f1.result()           
        df[['age','group_age']] = f2.result()          
        df['group_status'] =  f3.result()              
        df[['state_key','region_key']] = f4.result()  
        return df
    
    def __linage_normalized__(self, df:pd.DataFrame, dic:dict) -> list:
        variant = list()
        __pass = len(variant)
        for _df in df:
            if re.search('AY',_df):
                variant.append('Delta')
            else:
                for k,v in dic.items():
                    if _df in v:
                        variant.append(k)
                        break
                    if k == 'Other linages':
                        variant.append(k)
                        break
        if __pass == len(variant):
            logging.error("No se identifica la clasificación de: "+ _df)
        return variant
    
    def __get_quinquenios__(self,val):
        for k,v in age_quinquennia.items():
            if val in v:
                return (val,k)
    
    def __age_normalized__(self,df,dic)-> list:
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
            try:
               val = int(df['Patient age'].iloc[i])
               variant.append(self.__get_quinquenios__(val)) 
            except:
                if df['Patient age'].iloc[i].lower() in dic.keys():
                    val = dic[df['Patient age'].iloc[i].lower()]
                    variant.append(self.__get_quinquenios__(val))
                elif df['Patient age'].iloc[i] in ['Hospitalized', 'Ambulatory','Male','Female']:
                    variant.append(self.__get_quinquenios__(int(df['Gender'].iloc[i])))
                elif bool(re.search('A',df['Patient age'].iloc[i])):
                    val = re.search('A',df['Patient age'].iloc[i])
                    variant.append(self.__get_quinquenios__(int(df['Patient age'].iloc[i][:val.span()[0]])))
                elif bool(re.search('onths',df['Patient age'].iloc[i])):# referencia a meses
                    variant.append(self.__get_quinquenios__(0))
                elif bool(re.search('\.',df['Patient age'].iloc[i])):
                    val = re.search('\.',df['Patient age'].iloc[i])
                    variant.append(self.__get_quinquenios__(int(df['Patient age'].iloc[i][:val.span()[0]])))
                else:
                    variant.append(self.__get_quinquenios__(-1))
                    logging.error("No se identificó la edad:" + df['Patient age'].iloc[i])
        return variant
    
   

      
    def __state_normalized__(self,df,dic,dic_k) -> list:
        """_summary_
            clave y valor real del estado unificado
        Args:
            df (DataFrame): _description_
            dic (dict): _description_
        Returns:
            list<tuple>: clave y valor de estados 
        """   
        evaluate = lambda _df: (_df,dic_k[str(_df)]) if _df in dic_k.keys() else (99,"Extra")
        return [evaluate(dic[_df]) for _df in df] 

    
    def __state_group_patient_status__(self,df:pd.DataFrame, dic:dict) -> list:
        variant = list()
        for _df in df:
            try:
                variant.append(int(dic[_df.lower()]))
            except:
                if re.search('home',_df.lower()):
                    variant.append(1)
                else:
                    logging.error("No se identifica el estado:" + _df)
        return variant
    
    