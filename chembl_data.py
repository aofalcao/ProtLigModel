#from chembl_webresource_client.new_client import new_client
import pandas as pd
import requests, sys, json
import numpy as np
import math

WEBSITE_API = "https://rest.uniprot.org"

# Helper function to download data
def get_url(url, **kwargs):
    response = requests.get(url, **kwargs);

    if not response.ok:
        print(response.text)
        response.raise_for_status()
    sys.exit()

    return response


def get_uniprot_chembl_ids(uniprot_list):
    uniprot_chembl_dict = dict()
    chembl_ids_list = list()

    for prot in uniprot_list:
        r = get_url("https://www.uniprot.org/uniprotkb/{}.txt".format(prot))
        uniprot_res = r.text.splitlines()

        for l in uniprot_res:
            if 'ChEMBL' in l:
                chembl_ids = [i.strip() for i in l.split(';')[1:]]
                chembl_ids = [i for i in chembl_ids if '-.'not in i]
                #print(l_split)
                uniprot_chembl_dict[prot] = chembl_ids
                for i in chembl_ids:
                    chembl_ids_list.append(i)
        #print("https://www.uniprot.org/uniprotkb/{}.txt".format(l))

    return uniprot_chembl_dict, chembl_ids_list


file = open('uniprot-compressed_true_download_true_format_list_query__28gpcr_29_2-2022.10.06-10.13.04.74.list')
lines = file.readlines()
file.close()
uniprot_list = [l.strip() for l in lines]

ids_dict, chembl_list = get_uniprot_chembl_ids(uniprot_list)
chembl_list = list(set(chembl_list))




def read_url(url):
    try:
        response = requests.get(url)
        data = json.loads(response.text)
        return data
    except:
        return None
    


def read_activities(chembl_id, chembl_site = "https://www.ebi.ac.uk"):
    afields=['assay_type','assay_chembl_id','canonical_smiles','molecule_chembl_id','relation','standard_type',
    'target_organism','target_tax_id','type','units','value','standard_units','standard_upper_value',
                  'standard_value','standard_relation']
         #reads the activities of a given chembl_id (chembl target) url=
    url = chembl_site+"/chembl/api/data/activity.json?target_chembl_id="+chembl_id+"&limit=0"
    data=read_url(url)
    if data is not None:
        acts=data["activities"]

             ###############################################
             #removing these next lines will get at most 1000 molecules per problem. This is ok for testing ONLY
             # in production take the comments out
             ############################################
        while (data is not None) and (data["page_meta"] is not None) and (data["page_meta"]["next"] is not None):
            url = chembl_site+data["page_meta"]["next"]
            data=read_url(url)
            try:
                acts+=data["activities"]
             #print(url, len(acts))
            except:
                continue

        activities=[]

        for act in acts:
            #this Dictionary has been used in CQSAR and is therefore used here
            D={"assay_id": act["assay_chembl_id"], "target_id": act['target_chembl_id'],
               "assay_type": act['standard_type'], "units": act['standard_units'],
               "value": act['standard_value'], "rel":act['standard_relation'],
               "pvalue": act['pchembl_value'], "canonical_smiles": act['canonical_smiles'],
               "molecule_chembl_id": act['molecule_chembl_id'], "activity_comment": act['activity_comment'],
               "year": act['document_year'], "assay_description": act['assay_description']}

            activities.append(D)

        return activities
    else:
        print(data)
        return []
    



def get_all_activities(chembl_list):

    options = ['IC50', 'Ki']

    activities_list = list()
    without_activities = list()
    
    for target in chembl_list:
        res = read_activities(target, chembl_site = "https://www.ebi.ac.uk")
        if len(res)>0:
            target_activity_list = list()
            for r in res:
                if r['assay_type'] in options:
                    target_activity_list.append(r)
            if len(target_activity_list)>0:
                activities_list.append(target_activity_list)
            else:
                without_activities.append(target)
        else:
            without_activities.append(target)
    
    return activities_list, without_activities


activities, without_activities = get_all_activities(chembl_list)

ic50_list = []
ki_list = []

for acts in activities:
    df = pd.DataFrame(acts)
    df = df.astype({'value':'float', 'pvalue':'float', 'year':'float'})
    df_ki = df[df['assay_type']=='Ki'].copy()
    df_ic50 = df[df['assay_type']=='IC50'].copy()
    if len(df_ki)>0:
        ki_list.append(df_ki)
    if len(df_ic50)>0:
        ic50_list.append(df_ic50)







comments_list = ['Not Determined',
                 'Not Active (inhibition < 50% @ 10 uM and thus dose-reponse curve not measured)',
                 'ND(Insoluble)',
                 'Not Active',
                 'Nd(Insoluble)',
                 'No data',
                 'Non-toxic',
                 'NT',
                 'NC',
                 'Not Tested',
                 'Dose-dependent effect',
                 'Partial antagonist',
                 'Not Evaluated',
                 'No displacement',
                 'No effect']



def write_dict(mol_id, target_id, smiles, pvalue, obs):
    d = dict()
    d['molecule_chembl_id']=mol_id
    d['target_id']=target_id
    d['pvalue']=pvalue
    d['canonical_smiles']=smiles
    d['obs']=obs
    return d


def get_smaller_rel_obs(df_smaller): 
    arr = df_smaller.to_numpy()
    try:
        index = np.argmin(arr, axis=0)[4]
    except:
        index = 0
    obs = arr[index][5] + ' ' + str(arr[index][4])
    return obs


def one_mol(df_mol):
    d = dict()
    n = df_mol.to_numpy()
    obs=None
    pvalue=None
    if (n[0][9] is not None) and (n[0][9] in comments_list):
        obs = n[0][9]
        pvalue=0
    elif (n[0][9] is not None) and n[0][9] == 'Active':
        if np.isnan(n[0][6]):
            pvalue = 9 - math.log10(n[0][4])
        else:
            pvalue = n[0][6]
    elif (n[0][5]=='>') or (n[0][5]=='>='):
        pvalue=0
        obs=n[0][5] + ' ' + str(n[0][4])
    elif (n[0][5]=='<') or (n[0][5]=='<='):
        if np.isnan(n[0][6]):
            pvalue = 9 - math.log10(n[0][4])
            obs = n[0][5] + ' ' + str(n[0][4])
        else:
            pvalue=n[0][6]
            obs = n[0][5] + ' ' + str(n[0][4])
    elif not np.isnan(n[0][6]):
        pvalue=n[0][6]
    elif n[0][5]=='=' and (not np.isnan(n[0][4])):
        pvalue = 9 - math.log10(n[0][4])
      
    if pvalue is not None:
        d = write_dict(n[0][8], n[0][1], n[0][7], pvalue, obs)
    return d




def two_mol(df_mol):
    n = df_mol.to_numpy()
    d = dict()
    obs = None
    pvalue=None
    # se ambos tiverem pvalue
    if not(np.isnan(n[0][6])) and not(np.isnan(n[1][6])):
        # verificar se são muito diferentes
        if abs(n[0][6] - n[1][6]) < 1:
            pvalue = np.mean(n[:,6])
        # escolher o mais recente
        elif n[0][10] > n[1][10]:
            pvalue = n[0][6]
        else:
            pvalue = n[1][6]
        d = write_dict(n[0][8], n[0][1], n[0][7], pvalue, obs)
    # escolher o que tiver pvalue - significa que está ativo e rel é =
    # deixar na obs se houver rel < ou <=, a rever
    elif not(np.isnan(n[0][6])) or not(np.isnan(n[1][6])):
        #verificar os que são active mas não têm pvalue
        if ((n[0][9]=='Active' and np.isnan(n[0][6])) or (n[0][5]=='=' and (not np.isnan(n[0][4])))) and not(np.isnan(n[1][6])):
            active_pvalue = 9 - math.log10(n[0][4])
            if abs(active_pvalue - n[1][6]) < 1:
                pvalue = (active_pvalue + n[1][6])/2
            #em caso de diferença escolher o que já tinha pvalue
            elif n[0][10] > n[1][10]:
                pvalue = active_pvalue
            else:
                pvalue = n[1][6]

        elif ((n[1][9]=='Active' and np.isnan(n[1][6])) or (n[1][5]=='=' and (not np.isnan(n[1][4])))) and not(np.isnan(n[0][6])):
            active_pvalue = 9 - math.log10(n[1][4])
            if abs(active_pvalue - n[0][6]) < 1:
                pvalue = (active_pvalue + n[0][6])/2
            #em caso de diferença escolher o que já tinha pvalue
            elif n[0][10] < n[1][10]:
                pvalue = active_pvalue
            else:
                pvalue = n[0][6]

        # no caso de apenas um ter pvalue mas outro ter < ou <= na rel, a rever
        # de momento a guardaro pvalue do que tem pvalue e guardar o outro em obs
        elif not(np.isnan(n[0][6])):
            pvalue = n[0][6]
            if n[1][5]=='<' or n[1][5]=='<=':
                obs = n[1][5] + ' ' + str(n[1][4])

        else:
            pvalue = n[1][6]
            if n[0][5]=='<' or n[0][5]=='<=':
                obs = n[0][5] + ' ' + str(n[0][4])

        d = write_dict(n[0][8], n[0][1], n[0][7], pvalue, obs)
    
    elif n[0][9]=='Active' and np.isnan(n[0][6]) and n[1][9]=='Active' and np.isnan(n[1][6]):
        pvalue = (9 - math.log10(n[0][4]) + 9 - math.log10(n[1][4]))/2
        d = write_dict(n[0][8], n[0][1], n[0][7], pvalue, obs)

    elif n[0][9]=='Active' or n[1][9]=='Active':
        if n[0][9]=='Active':
            pvalue = 9 - math.log10(n[0][4])
        else:
            pvalue = 9 - math.log10(n[1][4])

        d = write_dict(n[0][8], n[0][1], n[0][7], pvalue, obs)

    else:
        if (n[0][5]=='<' or n[0][5]=='<=') and (n[1][5]=='<' or n[1][5]=='<='):
            if n[1][4] < n[0][4]:
                pvalue = 9 - math.log10(n[1][4])
                obs = n[1][5] + ' ' + str(n[1][4])
            else: 
                pvalue = 9 - math.log10(n[0][4])
                obs = n[0][5] + ' ' + str(n[0][4])
        elif (n[0][5]=='<' or n[0][5]=='<=') or (n[1][5]=='<' or n[1][5]=='<='):
            if n[0][5]=='<' or n[0][5]=='<=':
                pvalue = 9 - math.log10(n[0][4])
                obs = n[0][5] + ' ' + str(n[0][4])
            else:
                pvalue = 9 - math.log10(n[1][4])
                obs = n[1][5] + ' ' + str(n[1][4])
        else:
            if n[0][10] < n[1][10]:
                if n[1][9] is not None:
                    if n[1][9] in comments_list:
                        pvalue = 0
                        obs = n[1][9]
                else:
                    pvalue = 0
                    obs = n[1][5] + ' ' + str(n[1][4])
            else:
                if n[0][9] is not None:
                    if n[0][9] in comments_list:
                        pvalue = 0
                        obs = n[0][9]
                else:
                    pvalue = 0
                    obs = n[0][5] + ' ' + str(n[0][4])

        d = write_dict(n[0][8], n[0][1], n[0][7], pvalue, obs)

    return d




def curation(targets_activities):

    curated_list = []
    for i1,df in enumerate(targets_activities):
        #print(i1)
        p_test=False
        mols = list(df['molecule_chembl_id'].unique())
        if i1==204:
            print()
        l = []
        df = df.loc[~((df['activity_comment']=='Active') & (df['value'].isna()))
                   & ~((df['rel']=='=') & (df['value']<=0))
                   & ((~df['rel'].isna()) | (~df['activity_comment'].isna()))
                   & ((df['units']=='µM') | (df['units']=='nM'))]
        for i2,mol in enumerate(mols):
            df_mol = df.loc[(df['molecule_chembl_id']==mol)]  
                            #((~df['rel'].isna()) | (~df['activity_comment'].isna())) & 
                            #((df['units']=='µM') | (df['units']=='nM'))]
            
            #acertar unidades
            if len(df_mol[df_mol['units']=='µM']):
                cols = list(df_mol.columns)
                arr = df_mol.to_numpy()
                for index in range(arr.shape[0]):
                    if arr[index][3]=='µM':
                        arr[index][3] = 'nM'
                        arr[index][4] = arr[index][4] * 1000
                df_mol = pd.DataFrame(arr, columns=cols)
            if len(df_mol) == 1:
                d = one_mol(df_mol)
                if len(d)>0:
                    l.append(d)

            elif len(df_mol) == 2:
                d = two_mol(df_mol)
                if len(d)>0:
                    l.append(d)
            else:
                d = dict()
                pvalue_not_nan = df_mol[~np.isnan(df_mol['pvalue'])]
                active_pvalue_nan = df_mol[((np.isnan(df_mol['pvalue'])) & (df_mol['activity_comment']=='Active') & 
                                           (~np.isnan(df_mol['value']))) |
                                           ((np.isnan(df_mol['pvalue'])) & (df_mol['rel']=='=') & 
                                           (~np.isnan(df_mol['value'])))]
                smaller_rel_df = df_mol[(df_mol['rel'].str.contains('<')==True)]
                valid_len = len(pvalue_not_nan) + len(active_pvalue_nan)
                smaller_len = len(smaller_rel_df)
                pvalue = None
                obs = None
                # prioridade aos que têm pvalue ou são active
                if valid_len > 2:
                    active_pvalue_nan_np = None
                    active_bol = False
                    pvalue_not_nan_np = None
                    not_nan_bol = False
                    if len(active_pvalue_nan)>0:
                        active_pvalue_nan_np = active_pvalue_nan.to_numpy()
                        active_bol = True
                        for index in range(active_pvalue_nan_np.shape[0]):
                            active_pvalue_nan_np[index][6] = 9 - math.log10(active_pvalue_nan_np[index][4])
                    if len(pvalue_not_nan)>0:
                        pvalue_not_nan_np = pvalue_not_nan.to_numpy()
                        not_nan_bol = True
                    if not_nan_bol==True and active_bol==True:
                        n = np.concatenate((active_pvalue_nan_np,  pvalue_not_nan_np), axis=0)
                        pvalue = np.mean(n[:,6])
                    elif active_bol==False:
                        pvalue = np.mean(pvalue_not_nan_np[:,6])
                    else:
                        pvalue = np.mean(active_pvalue_nan_np[:,6])
                    if smaller_len>0:
                        obs = get_smaller_rel_obs(smaller_rel_df)
                    if active_bol:
                        n = active_pvalue_nan.to_numpy()
                        d = write_dict(n[0][8], n[0][1], n[0][7], pvalue, obs)
                    else:
                        n = pvalue_not_nan.to_numpy()
                        d = write_dict(n[0][8], n[0][1], n[0][7], pvalue, obs)

                elif valid_len == 2:
                    if len(pvalue_not_nan)==2:
                        d = two_mol(pvalue_not_nan)
                    elif len(active_pvalue_nan)==2:
                        d = two_mol(active_pvalue_nan)
                    else:
                        d = two_mol(pd.concat([pvalue_not_nan, active_pvalue_nan]))
                    if smaller_len>0:
                        d['obs'] = get_smaller_rel_obs(smaller_rel_df)

                elif valid_len == 1:
                    if len(pvalue_not_nan) == 1:
                        d = one_mol(pvalue_not_nan)
                    else:
                        active_pvalue_nan_cols = list(active_pvalue_nan.columns)
                        active_pvalue_nan_np = active_pvalue_nan.to_numpy()
                        active_pvalue_nan_np[0][6] = 9 - math.log10(active_pvalue_nan_np[0][4])
                        d = one_mol(pd.DataFrame(active_pvalue_nan_np, columns=active_pvalue_nan_cols))
                    if smaller_len>0:
                        d['obs'] = get_smaller_rel_obs(smaller_rel_df)

                else:
                    # se houver alguma em que rel é < ou <=, usar essa
                    #rever: média, mais recente ou valor mais baixo
                    if smaller_len>0:
                        smaller_np = smaller_rel_df.to_numpy()
                        try:
                            index = np.argmin(smaller_np, axis=0)[4]
                        except:
                            index = 0
                        pvalue = 9-math.log10(smaller_np[index][4])
                        obs = smaller_np[index][5] + ' ' + str(smaller_np[index][4])
                        d = write_dict(smaller_np[0][8], smaller_np[0][1], 
                                       smaller_np[0][7], pvalue, obs)
                    elif len(df_mol)>0:
                        pvalue = 0
                        n = df_mol.to_numpy()
                        try:
                            index = np.argmax(n, axis=0)[10]
                        except:
                            index = 0
                        if n[index][9] is not None:
                            obs = n[index][9]
                            d = write_dict(n[index][8], n[index][1], n[index][7], 0, obs)
                        else:
                            obs = n[index][5] + ' ' + str(n[index][4])
                            d = write_dict(n[index][8], n[index][1], n[index][7], 0, obs)
                if len(d)>0:
                    l.append(d)
        curated_list.append(pd.DataFrame(l))
        p_test = False
        #clear_output()
    return curated_list



ki_list_curated = curation(ki_list)

for i in ki_list_curated:
    if len(i)>0:
        id = i['target_id'].unique()[0]
        file = 'activities/curated/Ki/' + str(id) + '.csv'
        i.to_csv(file, sep=',', encoding='utf-8', index=False)



def function(value):
    if value=='active':
        value = 'Active'
    elif value=='inactive':
        value = 'Inactive'
    return value

for df in ic50_list:
    df['activity_comment'] = df['activity_comment'].apply(lambda x: function(x))
ic50_list_curated = curation(ic50_list)

for index,i in enumerate(ic50_list_curated):
    if len(i)>0:
        id = i['target_id'].unique()[0]
        file = 'activities/curated/IC50/' + str(id) + '.csv'
        i.to_csv(file, sep=',', encoding='utf-8', index=False)