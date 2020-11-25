import numpy as np
from keras.utils import Sequence

def shannon_entropy(x, tol=0.):
    return -np.sum(np.where(x > tol, x * np.log(x), 0), axis=-1)
    # return np.sum([-x * np.log(x) if x > 0 else 0 for _ in x])

def get_sparsity(x, tol=0.):
    return np.sum(x <= tol, axis=-1)/x.shape[1]
    
class DataLoader(Sequence):
    def __init__(self, dataset, batch_size, config, y_columns=None,
                 index=None, processor=None, shuffle=True, **kwargs):
        self.dataset = dataset
        if index is not None:
            self.dataset = self.dataset[index]
            assert len(index) == len(self.dataset)
        
        self.config = config
        self.batch_size = batch_size
        self.y_columns = y_columns
        self.shuffle = shuffle
        self.processor = processor
        
        self.on_epoch_end()
    
    def __len__(self):
        return int(np.ceil(1.0 * len(self.dataset) / self.batch_size))

    def __getitem__(self, index):
        """ Generate one batch of data. """
        # Generate indexes of the batch
        s = index * self.batch_size
        e = s + self.batch_size
        dataset = self.dataset[s:e]
        if self.y_columns is not None:
            X = dataset[:, :self.y_columns]
            y = dataset[:, self.y_columns:]
        else:
            X, y = dataset, None
        if self.processor:
            X, y = self.processor(X, y)
        
        return X, y

    def on_epoch_end(self):
        """ Updates indexes after each epoch. """
        if self.shuffle:
            np.random.shuffle(self.dataset)