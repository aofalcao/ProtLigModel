import numpy as np
import pickle
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, explained_variance_score

ids=[]
ki_uniprot_id_list = pickle.load(open('ki_uniprot_id_list.pickle', 'rb'))


#ids treino e teste2
#ids_test = ['O00222', 'P08172', 'P21728', 'P25024', 'P25098', 'P34947', 'P41594', 'Q8NFJ6']
ids_test = pickle.load(open('ids_test16_ki_d5_ordered.pickle', 'rb'))
ids_train = ki_uniprot_id_list
for i in ids_test:
    if i in ids_train:
        ids_train.remove(i)
ids_test = np.unique(ids_test)
ids_train = np.array(ids_train)   
       
dict_id_mol_transformed = pickle.load(open(f'dict_mol_id_ki.pickle', 'rb'))
dict_id_fp_transformed = pickle.load(open(f'dict_prot_id_transformed_ki_d5.pickle', 'rb'))

X_train = []
Y_train = []
mol_ids = []
for i in ids_train:
    mol_tuples = pickle.load(open('ki_mol_fp/{}_with_id.pickle'.format(i),'rb'))
    for mol_t in mol_tuples:
        # fingerprint: mol fingerprints + protein fingerprints
        mol_ids.append(mol_t[0])
        mol_fp = dict_id_mol_transformed[mol_t[0]]
        fp = mol_fp + dict_id_fp_transformed[i]
        X_train.append(fp)
        # guardar activity value
        Y_train.append(mol_t[2])

X_train = np.array(X_train)
Y_train = np.array(Y_train)

X_Y_test = []
mol_ids = []
for index, i in enumerate(ids_test):
    mol_ids = []
    X_Y_test.append([i, [], []])

    mol_tuples = pickle.load(open('ki_mol_fp/{}_with_id.pickle'.format(i),'rb'))
        
    for mol_t in mol_tuples:
        # fingerprint: mol fingerprints + protein fingerprints
        mol_ids.append(mol_t[0])
        mol_fp = dict_id_mol_transformed[mol_t[0]]
        fp = mol_fp + dict_id_fp_transformed[i]
        X_Y_test[index][1].append(fp)
        # guardar activity value
        X_Y_test[index][2].append(mol_t[2])     

    X_Y_test[index][1] = np.array(X_Y_test[index][1])
    X_Y_test[index][2] = np.array(X_Y_test[index][2])          

results = {'id':ids_test.tolist(), 'RMSE':[0]*len(ids_test), 'exp_var':[0]*len(ids_test)}


for it in range(10):
    print(it)
    model = RandomForestRegressor(n_estimators=200, n_jobs=-1).fit(X_train, Y_train)

    
    numerator_rmse = 0
    numerator_expvar = 0
    denominator = 0
    for index,t in enumerate(X_Y_test):
        if it==0:
            results['id'][index] = t[0]
       
       
        Y_test_pred = model.predict(t[1])
        rmse_test = mean_squared_error(t[2], Y_test_pred, squared=False) 
        expvar_test = explained_variance_score(t[2], Y_test_pred)
        results['RMSE'][index] += rmse_test/10
        results['exp_var'][index] += expvar_test/10


pickle.dump(results, open('ki/individual_results_mol_16_ki_d5.pickle', 'wb'))