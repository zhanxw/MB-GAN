import os
import numpy as np
import pandas as pd
import pickle
from keras.utils import Sequence

TOL = 1e-4

def shannon_entropy(x, tol=0.):
    return -np.sum(np.where(x > tol, x * np.log(x), 0), axis=-1)


def get_sparsity(x, tol=TOL):
    return np.sum(x <= tol, axis=-1)/x.shape[1]


def load_sample_pickle_data(filename="./data/raw_data.pkl"):
    ## Load raw dataset
    raw_data = pickle.load(open(filename, 'rb'))
    dataset = raw_data.iloc[:,1:].values/100.
    labels = raw_data["group"].values
    taxa_list = raw_data.columns[1:]
    
    data_o_case = dataset[labels == 'case']
    data_o_ctrl = dataset[labels == 'ctrl']
    
    return data_o_case, data_o_ctrl, taxa_list


def load_data_from_csv_file(csv_file, seq_dep=100, **kwargs):
    data = pd.read_csv(csv_file, **kwargs)
    ## extract species information only
    names = [_ for _ in data.index if _.split('|')[-1].startswith('s__')]
    # data_1_s = data.loc[data.index.str.contains('\|s__')].transpose()
    data_s = data.loc[names, :].transpose()
    print(data_s.shape)
    print(describe(shannon_entropy(data_s.values/seq_dep)))
    
    return data_s.values/seq_dep, data_s.columns


def expand_phylo(taxa_list):
    """ expand taxa to higher order. 
        adj_matrix, taxa_indices = expand_phylo(colnames)
    """
    memo, Ntaxa, adj, = {}, len(taxa_list), []
    for i, taxa in enumerate(taxa_list):
        memo[taxa] = i
        trunc = taxa.split('|')
        for j in range(len(trunc)):
            p = '|'.join(trunc[:j+1])
            if p not in memo:
                memo[p] = Ntaxa
                Ntaxa += 1
            adj.append((memo[taxa], memo[p]))
    # return adj, memo
    return adj, dict((v, k) for k, v in memo.items())


def adjmatrix_to_dense(x, shape, val=1):
    mask = np.zeros(shape)
    x = np.array(x).transpose()
    mask[tuple(x)] = val
    return mask
