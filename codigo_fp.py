import numpy as np
import pandas as pd
import scipy
import scipy.spatial
import pickle

def extractData(data):

    aas = {'ALA': 'A',
            'ARG': 'R',
            'ASN': 'N',
            'ASP': 'D',
            'ASX': 'B',
            'CYS': 'C',
            'GLU': 'E',
            'GLN': 'Q',
            'GLX': 'Z',
            'GLY': 'G',
            'HIS': 'H',
            'ILE': 'I',
            'LEU': 'L',
            'LYS': 'K',
            'MET': 'M',
            'PHE': 'F',
            'PRO': 'P',
            'SER': 'S',
            'THR': 'T',
            'TRP': 'W',
            'TYR': 'Y',
            'VAL': 'V'}
    #creades distance matrix, atom coordinates and atoms -> amino acids
    atoms_aas=[]
    coords=[]
    for line in data:
        #fs=line.strip().split()
        ATOM = line[:4]
        if ATOM=="ATOM":
            #ATOM, aid, atm, aa, lixo, aaid, x,y,z, s1, s2, atm2=fs
            aid = line[6:11].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            aa = line[17:20].strip()
            aaid = line[22:26].strip()
            aid=int(aid)
            #x,y,z, = float(x), float(y), float(z)
            atoms_aas.append(aas[aa]+aaid)
            coords.append([x,y,z])

    coords=np.array(coords)
    atoms_aas=np.array(atoms_aas)        
    DM = scipy.spatial.distance_matrix(coords, coords)
    return DM, coords, atoms_aas

def process_fp(fp):
    #returns processed aminoacids without ids from aas, but distinguishing different occurrences of the same AA
    L=sorted(list(fp))
    idx=0
    prev=L[0][0]
    L2=[prev+str(idx)]
    for i in L[1:]: 
        if i[0]==prev: idx+=1
        else: idx=0
        prev=i[0]
        L2.append(prev+str(idx))
    return frozenset(L2)

def compute_fps(DM, coords, atoms_aas, max_dist=5):
    fps=[]
    N=coords.shape[0]
    for i in range(N) :
        fp=set(atoms_aas[DM[i]<max_dist])
        if len(fp)>2:  
            fps.append(process_fp(fp))
    return set(fps)

def getFPs(fname, max_dist=5):
    data=open(fname, "rt").readlines()
    dm, crds, aaas = extractData(data)
    fps=compute_fps(dm, crds, aaas, max_dist)
    
    return fps

def h(fp_set, n=16381):
    fp = [0]*n
    hs = [hash(h) % n for h in fp_set]
    hs = set(hs)
    for item in hs:
        fp[item] = 1
    return fp

# get protein ids retrieved from uniprot
ids = []

with open('ids_144.txt', 'r') as fp:
    for line in fp:
        x = line[:-1]
        ids.append(x)

for i in ids:
    dist = 5
    fname="pdb_files/{}.pdb".format(i)
    fps= getFPs(fname, dist)
    fps = h(fps)
    pickle.dump(fps, open("fp_dist{}_ki/{}_D{}.pickle".format(dist,i,dist), "wb"))