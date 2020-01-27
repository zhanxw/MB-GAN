## Train a MB-GAN based on a sample case stool dataset
from model import MBGAN
from scipy.stats import describe
from utils import *

NAME = "mbgan_case"
EXP_DIR = "NielsenHB_2014_stool"
FILE = "./data/raw_data.pkl"
ITERS = 500000
BATCH_SIZE = 32
SAVE_INTERVAL = 1000


def get_save_fn(taxa_list):
    def fn(model, epoch):
        table_dir = os.path.join(model.log_dir, "tables")
        if not os.path.exists(table_dir):
            os.makedirs(table_dir)
        
        res = model.predict(1000, transform=None, seed=None)
        sparsity, entropy = get_sparsity(res), shannon_entropy(res)
        print("sparsity: %s" % str(describe(sparsity)))
        print("entropy: %s" % str(describe(entropy)))
        
        filename = "{:s}_{:06d}--{:.4f}--{:.4f}.csv".format(
            model.model_name, epoch, np.mean(sparsity), np.mean(entropy))
        
        pd.DataFrame(res, columns=taxa_list).to_csv(os.path.join(table_dir, filename))
        
    return fn


if __name__ == '__main__':
    # Load the dataset
    data_o_case, data_o_ctrl, taxa_list = load_sample_pickle_data(FILE)
    adj_matrix, taxa_indices = expand_phylo(taxa_list)
    tf_matrix = adjmatrix_to_dense(adj_matrix, shape=(len(taxa_list), len(taxa_indices)))
    
    model_config = {
        'ntaxa': 719,
        'latent_dim': 100,
        'generator': {'n_channels': 512},
        'critic': {'n_channels': 256, 'dropout_rate': 0.25, 
                   'tf_matrix': tf_matrix, 't_pow': 1000.}
    }
    
    train_config = {
        'generator': {'optimizer': ('RMSprop', {}), 'lr': 0.00005},
        'critic': {'loss_weights': [1, 1, 10], 
                   'optimizer': ('RMSprop', {}), 'lr': 0.00005},
    }
    
    ## Train for case:
    mbgan = MBGAN(NAME, model_config, train_config)
    mbgan.train(data_o_case, iteration=ITERS, batch_size=BATCH_SIZE, 
                n_critic=5, n_generator=1, save_fn=get_save_fn(taxa_list), 
                save_interval=SAVE_INTERVAL, experiment_dir=EXP_DIR,
                pre_processor=None, verbose=0)
    
