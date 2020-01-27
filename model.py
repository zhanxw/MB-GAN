# Large amount of credit goes to:
# https://github.com/keras-team/keras-contrib/blob/master/examples/improved_wgan.py

import keras.optimizers
import keras.backend as K

from keras.layers import Input, Dense, Dropout, Lambda, Layer
from keras.layers import BatchNormalization, Activation, LeakyReLU
from keras.layers.merge import _Merge
from keras.models import Sequential, Model

import os
import datetime
import numpy as np

from functools import partial


def gradient_penalty_loss(y_true, y_pred, averaged_samples):
    """ Calculates the gradient penalty loss for a batch of "averaged" samples. """
    gradients = K.gradients(y_pred, averaged_samples)[0]
    # compute the euclidean norm by squaring ...
    gradients_sqr = K.square(gradients)
    #   ... summing over the rows ...
    gradients_sqr_sum = K.sum(gradients_sqr,
                              axis=np.arange(1, len(gradients_sqr.shape)))
    #   ... and sqrt
    gradient_l2_norm = K.sqrt(gradients_sqr_sum)
    # compute lambda * (1 - ||grad||)^2 still for each single sample
    gradient_penalty = K.square(1 - gradient_l2_norm)
    # return the mean as loss over all the batch samples
    return K.mean(gradient_penalty)


def wasserstein_loss(y_true, y_pred):
    """ Calculates the Wasserstein loss for a sample batch. """
    return K.mean(y_true * y_pred)


def get_optimizer(optimizer, lr, decay=0.0, clipnorm=0.0, clipvalue=0.0, **kwargs):
    """ Get optimizer from keras.optimizers. """
    support_optimizers = {'SGD', 'RMSprop', 'Adagrad', 'Adadelta', 'Adam'}
    assert optimizer in support_optimizers
    fn = getattr(keras.optimizers, optimizer)
    return fn(lr, decay=decay, clipnorm=clipnorm, clipvalue=clipvalue, **kwargs)


class RandomWeightedAverage(_Merge):
    """ Calculate a random weighted average between two tensors. """
    def _merge_function(self, inputs):
        batch_size = K.shape(inputs[0])[0]
        alpha = K.random_uniform((batch_size, 1, 1, 1))
        return (alpha * inputs[0]) + ((1 - alpha) * inputs[1])


class PhyloTransform(Layer):
    def __init__(self, tf_matrix=None, **kwargs):
        if tf_matrix is None:
            self.kernel = None
        else:
            self.output_dim = tf_matrix.shape[1:]
            self.kernel = K.constant(tf_matrix, dtype='float32')
        super(PhyloTransform, self).__init__(**kwargs)

    def call(self, x):
        if self.kernel is None:
            return x
        else:
            return K.dot(x, self.kernel)
    
    def compute_output_shape(self, input_shape):
        if self.kernel is None:
            return input_shape
        else:
            return (input_shape[0], ) + self.output_dim


def build_generator(input_shape, output_units, n_channels=512):
    """ build the generator model. """
    model = Sequential()

    model.add(Dense(n_channels, activation="relu", input_shape=input_shape))
    model.add(BatchNormalization(momentum=0.8))
    model.add(Activation("relu"))
    model.add(Dense(n_channels))
    model.add(BatchNormalization(momentum=0.8))
    model.add(Activation("relu"))
    model.add(Dense(n_channels))
    model.add(BatchNormalization(momentum=0.8))
    model.add(Activation("relu"))
    model.add(Dense(output_units))
    model.add(Activation("softmax"))

    noise = Input(shape=input_shape)
    output = model(noise)

    return Model(noise, output)


def build_critic(input_shape, n_channels=256, dropout_rate=0.25, tf_matrix=None, t_pow=1000.):
    """ build the critic model. """
    model = Sequential()
    
    model.add(PhyloTransform(tf_matrix, input_shape=input_shape))
    model.add(Lambda(lambda x: K.log(1 + x * t_pow)/K.log(1 + t_pow)))
    model.add(Dense(n_channels))
    model.add(LeakyReLU(alpha=0.2))
    model.add(Dropout(dropout_rate))
    model.add(Dense(n_channels))
    model.add(LeakyReLU(alpha=0.2))
    model.add(Dropout(dropout_rate))
    model.add(Dense(n_channels))
    model.add(LeakyReLU(alpha=0.2))
    model.add(Dropout(dropout_rate))
    model.add(Dense(1))
    
    inputs = Input(shape=input_shape)
    validity = model(inputs)
    
    return Model(inputs, validity)


class MBGAN(object):
    def __init__(self, name, model_config, train_config):
        """ MBGAN model class. 
            name: provide a name for the given model/experiments.
            model_config: provide the configuration to build MBGAN:
                ntaxa: how many taxa are included in the real data
                latent_dim: the size of random vectors for the generator.
                generator: extra parameters parsed to build_generator.
                critic: extra parameters parsed to build_critic.
            train_config: provide the train configuration to build 
                computational graph. Includes: loss_weights, optimizer, 
                learning rate.
        """
        self.model_name = name
        self.model_config = model_config
        self.train_config = train_config
        
        self.ntaxa = self.model_config['ntaxa']
        self.latent_dim = self.model_config['latent_dim']
        
        # Build the generator and critic and construct the computational graph
        self.critic = build_critic((self.ntaxa,), 
                                   **self.model_config['critic'])
        self.generator = build_generator((self.latent_dim,), self.ntaxa, 
                                         **self.model_config['generator'])
        self.construct_critic_graph()
        self.construct_generator_graph()
    
    def construct_critic_graph(self):
        """ Construct computational graph for critic"""
        # Freeze generator's layers while training critic
        self.generator.trainable = False
        
        # Determines fake sample from given noise
        z = Input(shape=(self.latent_dim,))
        fake_sample = self.generator(z)
        fake = self.critic(fake_sample)
        
        # Determines real sample
        real_sample = Input(shape=(self.ntaxa,))
        valid = self.critic(real_sample)
        
        # Determines weighted average between real and fake sample
        interpolated_sample = RandomWeightedAverage()([real_sample, fake_sample])
        validity_interpolated = self.critic(interpolated_sample)
        
        # Get gradient penalty loss
        partial_gp_loss = partial(gradient_penalty_loss, averaged_samples=interpolated_sample)
        partial_gp_loss.__name__ = 'gradient_penalty'
        
        # Construct critic computational graph
        self.critic_graph = Model(inputs=[real_sample, z],
                                  outputs=[valid, fake, validity_interpolated])
        
        optimizer = get_optimizer(self.train_config['critic']['optimizer'][0], 
                                  lr=self.train_config['critic']['lr'], 
                                  **self.train_config['critic']['optimizer'][1])
        loss_weights = self.train_config['critic']['loss_weights']
        
        self.critic_graph.compile(
            loss=[wasserstein_loss, wasserstein_loss, partial_gp_loss],
            optimizer=optimizer, loss_weights=loss_weights,
        )
    
    
    def construct_generator_graph(self):
        """ Construct computational graph for generator. """
        # Freeze the critic's layers while training the generator
        self.critic.trainable = False
        self.generator.trainable = True
        
        # Generate sample and update generator
        z = Input(shape=(self.latent_dim,))
        fake_sample = self.generator(z)
        valid = self.critic(fake_sample)
        
        # Construct generator computational graph
        self.generator_graph = Model(z, valid)
        optimizer = get_optimizer(self.train_config['generator']['optimizer'][0], 
                                  lr=self.train_config['generator']['lr'], 
                                  **self.train_config['generator']['optimizer'][1])
        self.generator_graph.compile(loss=wasserstein_loss, optimizer=optimizer)
    

    def train(self, dataset, iteration, batch_size=32, n_critic=5, n_generator=1, 
              save_interval=50, save_fn=None, experiment_dir="mbgan_train",
              pre_processor=None, verbose=0, **kwargs):
        """ Train the MB-GAN with given dataset. 
            dataset: a sample x taxa tables.
            iteration: iteration to train the model.
            batch_size: samples to put in each batch.
            n_critic: times to update critic in each iteration.
            n_generator: times to update generator in each iteration.
            save_interval: the frequency to save model and sample outputs.
            save_fn: extra analyzing/saving functions run on save_interval.
            pre_processor: pre processing steps before feed into the model.
            post_processor: post processing steps on top of model generated results.
            verbose: not used. Aim to control screen output level.
            **kwargs: expanding
        """
        st_start = datetime.datetime.now()
        print("#####################################################")
        print("Training start at: {}".format(st_start.strftime('%Y-%m-%d %H:%M:%S')))
        print("Run MB-GAN for {:d} iterations with batch_size={:d}".format(iteration, batch_size))
        print("Save generated samples and model every {:d} iters".format(save_interval))
        
        # Create folders to save log file and models
        self.log_dir = os.path.join(experiment_dir, "{}_{:%Y%m%dT%H%M%S}".format(self.model_name, st_start))
        print("Results are exported to folder: {}".format(self.log_dir))
        if not os.path.exists(self.log_dir):
            print("    Create log folder: {}".format(self.log_dir))
            os.makedirs(self.log_dir)
        
        model_dir = os.path.join(self.log_dir, "models")
        if not os.path.exists(model_dir):
            print("    Create model folder: {}".format(model_dir))
            os.makedirs(model_dir)
        
        print("Generator structure:")
        self.generator.summary()
        print("Critic structure:")
        self.critic.summary()

        ## Set up the adversarial ground truths
        valid = -np.ones((batch_size, 1))
        fake =  np.ones((batch_size, 1))
        dummy = np.zeros((batch_size, 1)) # Dummy gt for gradient penalty
        
        ## Alternatively train critic (n_critic steps) and generator (n_genertor steps).
        for epoch in range(1, iteration+1):
            for _ in range(n_critic):
                # Randomly select a batch of samples to train the critic
                real = dataset[np.random.randint(0, dataset.shape[0], batch_size)]
                if pre_processor is not None:
                    real = pre_processor(real)
                noise = np.random.normal(0, 1, (batch_size, self.latent_dim))
                d_loss = self.critic_graph.train_on_batch([real, noise], [valid, fake, dummy])
            
            for _ in range(n_generator):
                #  Update the generator
                noise = np.random.normal(0, 1, (batch_size, self.latent_dim))
                g_loss = self.generator_graph.train_on_batch(noise, valid)

            # Plot the progress
            log_info = [
                "iter={:d}".format(epoch), 
                "[D loss={:.6f}, w_loss_real={:.6f}, w_loss_fake={:.6f}, gp_loss={:.6f}]".format(*d_loss),
                "[G loss={:.6f}]".format(g_loss),
            ]
            print("{} {} {}".format(*log_info))

            # Save generated samples on saving interval
            if epoch % save_interval == 0:
                print("At iter={:d}, saving model weights and exporting generated samples".format(epoch))
                # custom_objects={'PhyloTransform': PhyloTransform})
                self.critic.save(os.path.join(model_dir, "{:s}_{:06d}_critic.h5".format(self.model_name, epoch)))
                self.generator.save(os.path.join(model_dir, "{:s}_{:06d}_generator.h5".format(self.model_name, epoch)))
                save_fn(self, epoch)
    
    def predict(self, n_samples=100, transform=None, seed=None):
        np.random.seed(seed)
        z = np.random.normal(0, 1, (n_samples, self.latent_dim))
        res = self.generator.predict(z)
        if transform is not None:
            res = transform(res)
        
        return res
