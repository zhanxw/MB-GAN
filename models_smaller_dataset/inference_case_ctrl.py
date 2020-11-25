import os
import pickle
import argparse
import numpy as np
import pandas as pd
from scipy.stats import describe
from keras.models import load_model

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = ""

SEED = 256
TOL = 1e-4
GENERATOR_CASE_PATH = "./micro_zeller_case/micro_zeller_case_20191127T174114/models/micro_zeller_case_537000_generator.h5"
GENERATOR_CTRL_PATH = "./micro_zeller_ctrl/micro_zeller_ctrl_20191127T174509/models/micro_zeller_ctrl_559000_generator.h5"


def get_parser():
    ## Get command line parameters
    parser = argparse.ArgumentParser(description='Generate microbiome data')    
    parser.add_argument('--num_samples', default=1000, type=int, help='number of samples to generate')
    parser.add_argument('--generator_case_path', default=GENERATOR_CASE_PATH, type=str, help='generator case weights path')
    parser.add_argument('--generator_ctrl_path', default=GENERATOR_CTRL_PATH, type=str, help='generator ctrl weights path')
    parser.add_argument('--output_folder', default='outputs', type=str, help='generator weights path')
    parser.add_argument('--tol', default=TOL, type=float, help='ignore values smaller than this value')
    
    return parser


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

def load_data():
    data_1 = pd.read_csv('./data/MBGANBiomarker/zeller.ctrl.metaphlan_bugs_list.stool.tsv', sep='\t')
    data_2 = pd.read_csv('./data/MBGANBiomarker/zeller.case.metaphlan_bugs_list.stool.tsv', sep='\t')
    ## extract species information only
    name_1 = [_ for _ in data_1.index if _.split('|')[-1].startswith('s__')]
    data_1_s = data_1.loc[name_1, :].transpose()
    name_2 = [_ for _ in data_2.index if _.split('|')[-1].startswith('s__')]
    data_2_s = data_2.loc[name_2, :].transpose()
    
    data_0_s = pd.merge(data_1_s, data_2_s, 'outer')
    labels_0 = np.array(['ctrl'] * len(data_1_s) + ['case'] * len(data_2_s))
    # labels_0 = np.concatenate([np.zeros((len(data_1_s),1)), np.ones((len(data_2_s),1))])
    
    dataset, colnames = (data_0_s.values, labels_0), data_0_s.columns
    return dataset, colnames


if __name__ == "__main__":
    ## Load raw dataset
    opt = get_parser().parse_args()
    print(opt)
    (dataset, labels), taxa_list = load_data()
    data_o_case = dataset[labels == 'case']
    data_o_ctrl = dataset[labels == 'ctrl']

    ## Generate data
    print("Load case generator from: ")
    print(opt.generator_case_path)
    generator_case = load_model(opt.generator_case_path)
    
    print("Load ctrl generator from: ")
    print(opt.generator_ctrl_path)
    generator_ctrl = load_model(opt.generator_ctrl_path)
    
    print("Simulate %d case sample and %d ctrl sample" % (opt.num_samples, opt.num_samples))
    data_g_case = simulate(generator_case, n_samples=opt.num_samples, seed=SEED)
    data_g_ctrl = simulate(generator_ctrl, n_samples=opt.num_samples, seed=SEED)
    
    ## Show data statistics
    print("Sparsity")
    print(pd.DataFrame(
        [describe(get_sparsity(data_o_ctrl, opt.tol)),
         describe(get_sparsity(data_g_ctrl, opt.tol)),
         describe(get_sparsity(data_o_case, opt.tol)),
         describe(get_sparsity(data_g_case, opt.tol)),], 
        index=['Original ctrl', 'GAN ctrl', 'Original case', 'GAN case']))
    print("Shannon Entropy")
    print(pd.DataFrame(
        [describe(shannon_entropy(data_o_ctrl)),
         describe(shannon_entropy(data_g_ctrl)),
         describe(shannon_entropy(data_o_case)),
         describe(shannon_entropy(data_g_case)),], 
        index=['Original ctrl', 'GAN ctrl', 'Original case', 'GAN case']))
    
    
    ## Save simulated species level data
    if not os.path.exists(opt.output_folder):
        os.makedirs(opt.output_folder)
        
    pd.DataFrame(data_g_case, columns=taxa_list).to_csv(
        os.path.join(opt.output_folder, "zeller.case.metaphlan_bugs_list.stool.simu.taxa.csv"))
    pd.DataFrame(data_g_ctrl, columns=taxa_list).to_csv(
        os.path.join(opt.output_folder, "zeller.ctrl.metaphlan_bugs_list.stool.simu.taxa.csv"))
    
    ## Save simulated phylogenetic data
    adj_matrix, taxa_indices = expand_phylo(taxa_list)
    tf_matrix = adjmatrix_to_dense(adj_matrix, shape=(len(taxa_list), len(taxa_indices)))
    
    pd.DataFrame(np.dot(data_g_case, tf_matrix), columns=taxa_indices).to_csv(
        os.path.join(opt.output_folder, "zeller.case.metaphlan_bugs_list.stool.simu.phylo.csv"))
    pd.DataFrame(np.dot(data_g_ctrl, tf_matrix), columns=taxa_indices).to_csv(
        os.path.join(opt.output_folder, "zeller.ctrl.metaphlan_bugs_list.stool.simu.phylo.csv"))
