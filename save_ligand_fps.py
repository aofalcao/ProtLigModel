import numpy as np
import pandas as pd
import pickle
import requests, sys, json


# get pandas dataframe from filepath (activities curated)
def get_file(filepath):
    
    def function(value):
        if value<=5:
            value = 0
        elif value>=9:
            value = 1
        else:
            value = value/4 - 5/4
        return value
    
    file = pd.read_csv(filepath)
    file = file.astype({'pvalue':'float'})
    file['pvalue'] = file['pvalue'].apply(lambda x: function(x))
    file = file.drop(['target_id', 'obs'], axis=1)
    file = file.drop(file[file['canonical_smiles'].isna()].index)
    file = file.drop(file[file['pvalue'].isna()].index)
    
    return file



def save_fp_values(dict_uniprot_chembl, activities_directory='F:/Jupyter_directory/TESE/activities/curated/Ki/'):
    
    for uniprot_id in dict_uniprot_chembl:
        chembl_id = dict_uniprot_chembl[uniprot_id]
        # nem todos tem atividades
        try:
            file = get_file(activities_directory+'{}.csv'.format(chembl_id))
            file = file.to_numpy()
            # lista com tuplos: (fingerprint, activity_value)
            tuple_list = []

            for line in file:
                m = Chem.MolFromSmiles(line[2])
                fp = AllChem.GetMorganFingerprintAsBitVect(m,3,nBits=2048)
                fp = list(fp)
                tuple_fp_value = (line[0], fp, line[1])
                tuple_list.append(tuple_fp_value)

            pickle.dump(tuple_list, open("ki_mol_fp([fp],activity_value)/{}_with_id.pickle".format(uniprot_id), "wb"))
        

        except Exception as e:
            continue


dict_uniprot_chembl = pickle.load(open("dict_uniprot_chembl.pickle", "rb"))
save_fp_values(dict_uniprot_chembl)