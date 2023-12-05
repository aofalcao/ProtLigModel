import pickle
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error, explained_variance_score
import matplotlib.pyplot as plt
import seaborn as sns

dict_dimensions = pickle.load(open('dict_dimensions.pickle', 'rb'))

r = pickle.load(open('mol_global_results_dict_ki_general_d5_individual.pickle', 'rb'))
r_df = pd.DataFrame.from_dict(r)

ids = pickle.load(open('ki_uniprot_id_list.pickle', 'rb'))
#ids.remove('P21462')
#ids.remove('Q14833')

d = {'id':ids, 'RMSE_unified':[], 'RVE_unified':[], 'n_train_unified':[], 'n_test_unified':[]}
n_dict = {'id':[], 'n_train':[], 'n_test':[], 'run':[]}
   
    
    
for index,i in enumerate(ids):
    #if i==ids[index]:
    #    print(i)
    if i!=ids[index]:
        print('wrong: ', i)
        
    rmse = 0
    exp_var = 0
    total_n = 0
    
    for run in range(10):
        values = r_df[(r_df['run']==run) & (r_df['id']==i)][['real_value', 'pred_value']].to_numpy()
        real = values[:, 0]
        pred = values[:, 1]
        #n_dict['id'].append(i); n_dict['n'].append(len(real)); n_dict['run'].append(run)
        if len(values)>0:

            rmse += mean_squared_error(real, pred, squared=False) * len(real)
            exp_var += explained_variance_score(real, pred) * len(real)
            total_n += len(real)

            n_dict['id'].append(i)
            n_dict['n_train'].append(dict_dimensions[i]-len(real))
            n_dict['n_test'].append(len(real))
            n_dict['run'].append(run)

            n_train = dict_dimensions[i]-len(real)
            n_test = len(real)
            
        
    d['RMSE_unified'].append(rmse/total_n)
    d['RVE_unified'].append(exp_var/total_n)
    d['n_train_unified'].append(n_train)
    d['n_test_unified'].append(n_test)  
    
    
df_unified = pd.DataFrame.from_dict(d)