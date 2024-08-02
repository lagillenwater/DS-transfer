from utils import *
import pandas as pd

parser = define_arguments()
train_file1,test_file1, train_file2,test_file2,latent_dim, epochs, batch_size,output_dir = generate_arguments(parser)
vae1,encoder1,decoder1,history1 = transformer_vae(train_file1, test_file1,train_file2,test_file2, latent_dim,epochs, batch_size)

train_plot(history1,train_file1,output_dir)



