
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import keras
from keras import layers
from keras import backend as K
from utils import *


#Define arguments for each required and optional input
def define_arguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ## Required inputs
    parser.add_argument("--data-file1",dest="data_file1",required=True,help="full dataset 1")        
    parser.add_argument("--train-file1",dest="train_file1",required=True,help="training dataset 1")
    parser.add_argument("--test-file1",dest="test_file1",required=True,help="testing dataset 1")
    parser.add_argument("--output-dir",dest="output_dir",required=True,help="output directory")

    parser.add_argument("--data-file2",dest="data_file2",required=True,help="full dataset 2")        
    parser.add_argument("--train-file2",dest="train_file2",required=True,help="training dataset 2")
    parser.add_argument("--test-file2",dest="test_file2",required=True,help="testing dataset 2")

    # optional inputs
    parser.add_argument("--epochs", dest = "epochs", required=False, help="epochs", default=100, type = int)
    parser.add_argument("--latent-dim",dest="latent_dim",required=False,help="latent dimensions",default = 16, type = int)
    parser.add_argument("--batch-size",dest="batch_size",required=False,default=1000,help="batch size", type = int)
    return parser

# generate arguments
def generate_arguments(parser):

    #Generate argument parser and define arguments
    parser = define_arguments()
    args = parser.parse_args()

    data_file1 = args.train_file1
    train_file1 = args.train_file1
    test_file1 = args.test_file1

    data_file2 = args.train_file2
    train_file2 = args.train_file2
    test_file2 = args.test_file2    

    latent_dim = args.latent_dim
    epochs = args.epochs
    batch_size = args.batch_size
    output_dir = args.output_dir
    return data_file1,train_file1,test_file1,data_file2,train_file2,test_file2,latent_dim, epochs,batch_size,output_dir

# read in the training and testing files
def read_input_files( train_file,test_file):
    train = pd.read_csv(train_file, index_col = 0)
    test = pd.read_csv(test_file, index_col = 0)
    return train,test

# find the training shape
def get_input_dim(train):
    input_dim=train.shape[1]
    return input_dim



# output training plot
def train_plot(history,train_file,output_dir):
    plt.figure(figsize=(10, 5))
    plt.plot(history["loss"], label="Training data")
    plt.plot(history["val_loss"], label="Testing data")
    plt.ylabel("MSE + KL Divergence")
    plt.xlabel("No. Epoch")
    plt.legend()
    #plt.show()
    plt.savefig(output_dir +train_file.split("/")[-1].split(".")[0]+ 'training_plot.png')




def vanilla_vae(train_file, test_file, latent_dim, epochs, batch_size):
    x_train,x_test= read_input_files(train_file,test_file)

    original_dim = get_input_dim(x_train)
    intermediate_dim = 16
    latent_dim = 4

    inputs = keras.Input(shape=(original_dim,))
    h = layers.Dense(intermediate_dim, activation='relu')(inputs)
    z_mean = layers.Dense(latent_dim)(h)
    z_log_sigma = layers.Dense(latent_dim)(h)


    def sampling(args):
        z_mean, z_log_sigma = args
        epsilon = K.random_normal(shape=(K.shape(z_mean)[0], latent_dim),
                                  mean=0., stddev=0.1)
        return z_mean + K.exp(z_log_sigma) * epsilon

    z = layers.Lambda(sampling)([z_mean, z_log_sigma])

    # Create encoder
    encoder = keras.Model(inputs, [z_mean, z_log_sigma, z], name='encoder')

    # Create decoder
    latent_inputs = keras.Input(shape=(latent_dim,), name='z_sampling')
    x = layers.Dense(intermediate_dim, activation='relu')(latent_inputs)
    outputs = layers.Dense(original_dim, activation='sigmoid')(x)
    decoder = keras.Model(latent_inputs, outputs, name='decoder')

    # instantiate VAE model
    outputs = decoder(encoder(inputs)[2])
    vae = keras.Model(inputs, outputs, name='vae_mlp')

    reconstruction_loss = keras.losses.binary_crossentropy(inputs, outputs)
    reconstruction_loss *= original_dim
    kl_loss = 1 + z_log_sigma - K.square(z_mean) - K.exp(z_log_sigma)
    kl_loss = K.sum(kl_loss, axis=-1)
    kl_loss *= -0.5
    vae_loss = K.mean(reconstruction_loss + kl_loss)
    vae.add_loss(vae_loss)
    vae.compile(optimizer='adam')

    # x_train = x_train.astype('float32') / 255.
    # x_test = x_test.astype('float32') / 255.
    # x_train = x_train.reshape((len(x_train), np.prod(x_train.shape[1:])))
    # x_test = x_test.reshape((len(x_test), np.prod(x_test.shape[1:])))

    vae.fit(x_train,
            epochs=epochs,
            batch_size=batch_size,
            validation_data=(x_test,None))

    history = pd.DataFrame(vae.history.history)

    return vae,encoder,decoder,history
    
