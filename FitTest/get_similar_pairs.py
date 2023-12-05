import glob
import pickle
import numpy as np
import pandas as pd


ki_list = pickle.load(open('ki_uniprot_id_list.pickle', 'rb'))
total_activities_dict = pickle.load(open('total_activities_dict.pickle', 'rb'))

total_activities_dict_df = pd.DataFrame.from_dict(total_activities_dict)

index_dict = dict()
for index, i in enumerate(total_activities_dict['id']):
    if i in ki_list:
        index_dict[i] = index


filtered=total_activities_dict_df[total_activities_dict_df['ki']>=100]

ids_with_atleast_100 = list(filtered['id'].values)

all_data=[]
for i in ids_with_atleast_100:
    all_data.append(pickle.load(open(f'F:/Jupyter_directory/TESE/proteins_Secondary_structure/Código_e_dados/fp_dist5/{i}_D5.pickle', "rb")))


pairs = []
js=[]
for i, s1 in enumerate(all_data):
    for j, s2 in enumerate(all_data):
        jacc=len(s1 & s2)/len(s1 | s2)
        if jacc > .05 and i != j: 
            #print(i, j, "--->", jacc, len(s1 & s2), len(s1 | s2), ki_list[i], ki_list[j])
            if ((ids_with_atleast_100[i],ids_with_atleast_100[j]) not in pairs) and ((ids_with_atleast_100[j],ids_with_atleast_100[i]) not in pairs):
                pairs.append((ids_with_atleast_100[i], ids_with_atleast_100[j]))
                js.append(jacc)


js,pairs = zip(*sorted(zip(js,pairs)))

ids = list(pairs)
js2 = list(js)
ids_train = []
ids_test = []
tg = []
index = 0
for i in range(-1, -len(ids)-1, -1):
    # prot com mais atividades vai para treino
    index1 = index_dict[ids[i][0]]
    index2 = index_dict[ids[i][1]]
    
    if total_activities_dict['ki'][index1] >= total_activities_dict['ki'][index2]:
        if ids[i][1] not in ids_train and ids[i][0] not in ids_test:
            ids_train.append(ids[i][0])
            ids_test.append(ids[i][1])
            tg.append([ids[i][0], ids[i][1], js2[i]])
    else:
        if ids[i][0] not in ids_train and ids[i][1] not in ids_test:
            ids_train.append(ids[i][1])
            ids_test.append(ids[i][0])
            tg.append([ids[i][1], ids[i][0], js2[i]])
    if len(tg)==16:
        break
    
pickle.dump(ids_train, open('ids_train16_ki_d5_ordered.pickle','wb'))
pickle.dump(ids_test, open('ids_test16_ki_d5_ordered.pickle','wb'))





dict_ids_sim_n_activities = {'id1':[],'n_activities_ki_1':[],
                             'id2':[],'n_activities_ki_2':[],'similarity':[]}
for i in tg:
    dict_ids_sim_n_activities['id1'].append(i[0])
    dict_ids_sim_n_activities['id2'].append(i[1])
    dict_ids_sim_n_activities['similarity'].append(i[2])
    
    index1 = total_activities_dict['id'].index(i[0])
    index2 = total_activities_dict['id'].index(i[1])
    
    dict_ids_sim_n_activities['n_activities_ki_1'].append(total_activities_dict['ki'][index1])
    dict_ids_sim_n_activities['n_activities_ki_2'].append(total_activities_dict['ki'][index2])


print('average similarity: ', np.average(dict_ids_sim_n_activities['similarity']))
print('min: ', dict_ids_sim_n_activities['similarity'][-1])
print('max: ', dict_ids_sim_n_activities['similarity'][0])

pd.DataFrame.from_dict((dict_ids_sim_n_activities)).to_csv('new_pairs_d5.csv', index=False)