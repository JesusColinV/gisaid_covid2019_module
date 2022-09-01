from diasig.core.Diasing import Analyzer
from diasig.utils.utils import Counter
import pandas as pd

path = "gisaid_hcov-19_2022_08_17_04.tsv"

an = Analyzer(path=path)
df = an.load()
val = 0
#    val,columns = Counter.get_table(amino,)
#    for segment in amino.keys():
#        myDF = pd.DataFrame(val,columns=columns)
#    for i in range(len(val)-1):
#        myDF[i]['area'] = [x[0] for num in myDF[i].position for x in parts_protein[myDF[i].full[0]] if num >= x[1] and num <= x[2]]
#    myDF = pd.DataFrame()
#    
#    df_mutations_total = pd.concat([i.rename(columns = {'Spike_':'mutation','E_':'mutation','N_':'mutation','M_':'mutation','NSP':'mutation'}) for i in myDF]).set_index('mutation')
#    amino = Counter.get_mutations(df.query('Age >= 18')) 
#    myDF = Counter.get_table(amino)
#    df_mutations_over = Counter.pd.concat([i.rename(columns = {'Spike_':'mutation','E_':'mutation','N_':'mutation','M_':'mutation','NSP':'mutation'}) for i in myDF])
#    amino = Counter.get_mutations(df.query('Age < 18')) 
#    myDF = Counter.get_table(amino)
#    df_mutations_under = pd.concat([i.rename(columns = {'Spike_':'mutation','E_':'mutation','N_':'mutation','M_':'mutation','NSP':'mutation'}) for i in myDF])
#
#