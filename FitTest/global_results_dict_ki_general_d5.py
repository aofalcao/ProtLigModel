import numpy as np
import pickle
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, explained_variance_score


ki_uniprot_id_list = pickle.load(open('ki_uniprot_id_list.pickle', 'rb'))
    
    
results = {'RMSE':[], 'exp_var':[]}
individual_results = {'id':[], 'real_value':[], 'pred_value':[], 'run':[]}

dict_id_mol_transformed = pickle.load(open(f'dict_mol_id_ki.pickle', 'rb'))
dict_id_fp_transformed = pickle.load(open(f'dict_prot_id_transformed_ki_d5.pickle', 'rb'))
        
X = []
Y = []
ids = []
for i in ki_uniprot_id_list:
    mol_tuples = pickle.load(open('ki_mol_fp/{}_with_id.pickle'.format(i),'rb'))
    for mol_t in mol_tuples:
        # fingerprint: mol fingerprints + protein fingerprints
        ids.append(i)
        mol_fp = dict_id_mol_transformed[mol_t[0]]
        fp = mol_fp + dict_id_fp_transformed[i]
        X.append(fp)
        # guardar activity value
        Y.append(mol_t[2])
        
for run in range(10):
    print(run)
    X_train, X_test, Y_train, Y_test, ids_train, ids_test = train_test_split(X, Y, ids, test_size=0.2)
    X_train = np.array(X_train, dtype=np.float32)
    X_test = np.array(X_test, dtype=np.float32)
    Y_train = np.array(Y_train)
    Y_test = np.array(Y_test)

    model = RandomForestRegressor(n_estimators=200, n_jobs=-1).fit(X_train, Y_train)
    
    Y_test_pred = model.predict(X_test)
    rmse_test = mean_squared_error(Y_test, Y_test_pred, squared=False) 
    expvar_test = explained_variance_score(Y_test, Y_test_pred)
    
    for index in range(len(Y_test)):
        individual_results['id'].append(ids_test[index])
        individual_results['real_value'].append(Y_test[index])
        individual_results['pred_value'].append(Y_test_pred[index])
        individual_results['run'].append(run)

    print(rmse_test, expvar_test)
    results['RMSE'].append(rmse_test)
    results['exp_var'].append(expvar_test)
        
pickle.dump(results, open('ki/mol_global_results_dict_ki_general_d5.pickle', 'wb'))
pickle.dump(individual_results, open('ki/mol_global_results_dict_ki_general_d5_individual.pickle', 'wb'))