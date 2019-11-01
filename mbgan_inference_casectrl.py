import os
import pickle
import numpy as np
import pandas as pd
from scipy.stats import describe
from keras.models import load_model

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = ""

SEED = 256
TOL = 1e-4
GENERATOR_CASE_PATH = os.path.join('models', 'stool_2_case_generator.h5')
GENERATOR_CTRL_PATH = os.path.join('models', 'stool_2_ctrl_generator.h5')


def shannon_entropy(x, tol=0.):
    return -np.sum(np.where(x > tol, x * np.log(x), 0), axis=-1)


def get_sparsity(x, tol=0.):
    return np.sum(x <= tol, axis=-1)/x.shape[1]


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
    
    return adj, dict((v, k) for k, v in memo.items())


def adjmatrix_to_dense(x, shape, val=1):
    mask = np.zeros(shape)
    x = np.array(x).transpose()
    mask[tuple(x)] = val
    return mask


def simulate(model, n_samples=1000, transform=None, seed=None):
    np.random.seed(seed)
    latent_dim = model.inputs[0].shape[-1]
    z = np.random.normal(0, 1, (n_samples, latent_dim))
    res = model.predict(z)
    if transform is not None:
        res = transform(res)
    
    return res


if __name__ == "__main__":
    ## Load raw dataset
    raw_data = pickle.load(open("./data/raw_data.pkl", 'rb'))
    dataset = raw_data.iloc[:,1:].values/100.
    labels = raw_data["group"].values
    taxa_list = raw_data.columns[1:]
      
    data_o_case = dataset[labels == 'case']
    data_o_ctrl = dataset[labels == 'ctrl']

    ## Generate data
    generator_case = load_model(GENERATOR_CASE_PATH)
    generator_ctrl = load_model(GENERATOR_CTRL_PATH)
    data_g_case = simulate(generator_case, n_samples=1000, seed=SEED)
    data_g_ctrl = simulate(generator_ctrl, n_samples=1000, seed=SEED)

    ## Show data statistics
    print("Sparsity")
    print(pd.DataFrame(
        [describe(get_sparsity(data_o_ctrl, TOL)),
         describe(get_sparsity(data_g_ctrl, TOL)),
         describe(get_sparsity(data_o_case, TOL)),
         describe(get_sparsity(data_g_case, TOL)),], 
        index=['Original ctrl', 'GAN ctrl', 'Original case', 'GAN case']))
    print("Shannon Entropy")
    print(pd.DataFrame(
        [describe(shannon_entropy(data_o_ctrl)),
         describe(shannon_entropy(data_g_ctrl)),
         describe(shannon_entropy(data_o_case)),
         describe(shannon_entropy(data_g_case)),], 
        index=['Original ctrl', 'GAN ctrl', 'Original case', 'GAN case']))
    
    ## Save simulated species level data
    pd.DataFrame(data_g_case, columns=taxa_list).to_csv("./outputs/stools_2_case_species.csv")
    pd.DataFrame(data_g_ctrl, columns=taxa_list).to_csv("./outputs/stools_2_ctrl_species.csv")
    
    ## Save simulated phylogenetic data
    adj_matrix, taxa_indices = expand_phylo(taxa_list)
    tf_matrix = adjmatrix_to_dense(adj_matrix, shape=(len(taxa_list), len(taxa_indices)))
    
    pd.DataFrame(np.dot(data_g_case, tf_matrix), columns=taxa_indices).to_csv("./outputs/stools_2_case_phylo.csv")
    pd.DataFrame(np.dot(data_g_ctrl, tf_matrix), columns=taxa_indices).to_csv("./outputs/stools_2_ctrl_phylo.csv")
    
    
    
    
    