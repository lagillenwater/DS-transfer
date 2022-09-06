from vae import VAE
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np

def main():
#Define arguments for each required and optional input
    def define_arguments():
        parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        ## Required inputs
        parser.add_argument("--data-file",dest="data_file",required=True,help="full dataset")        
        parser.add_argument("--train-file",dest="train_file",required=True,help="training dataset")
        parser.add_argument("--test-file",dest="test_file",required=True,help="testing dataset")
        parser.add_argument("--output-dir",dest="output_dir",required=True,help="output directory")
        parser.add_argument("--encoder-architecture",dest="encoder_architecture",required=True,help="encoder architecture", nargs= '*', type = int, default = [])
        parser.add_argument("--decoder-architecture",dest="decoder_architecture",required=True,help="decoder architecture", nargs='*', type = int, default = [])

        # optional inputs
        parser.add_argument("--epochs", dest = "epochs", required=False, help="epochs", default=50, type = int)
        parser.add_argument("--latent-dim",dest="latent_dim",required=False,help="latent dimensions",default = 16, type = int)
        parser.add_argument("--batch-size",dest="batch_size",required=False,default=1000,help="batch size")
        parser.add_argument("--optimizer",dest="optimizer",required=False,default='adam',help="optimizer")
        parser.add_argument("--learning-rate",dest="learning_rate",required=False,default=0.0005,help="learning_rate", type=float)
        parser.add_argument("--epsilon",dest="epsilon_std",required=False,default=1.0,help="epsilon", type = float)
        parser.add_argument("--lambda",dest="lam",required=False,default=0.0,help="lambda", type = float)
        parser.add_argument("--beta",dest="beta",required=False,default=1.0,help="beta", type = float)
        parser.add_argument("--loss",dest="loss",required=False,default="binary_crossentropy",help="loss function")
        parser.add_argument("--encoder_batch_norm",dest="encoder_batch_norm",required=False,default=True,help="encoder batch norm")
        parser.add_argument("--verbose",dest="verbose",required=False,default=True,help="verbose")
        return parser

    # generate arguments
    def generate_arguments(parser):

        #Generate argument parser and define arguments
        parser = define_arguments()
        args = parser.parse_args()

        data_file = args.train_file
        train_file = args.train_file
        test_file = args.test_file
        latent_dim = args.latent_dim
        epochs = args.epochs
        batch_size = args.batch_size
        optimizer = args.optimizer
        learning_rate = args.learning_rate
        epsilon_std = args.epsilon_std
        beta =args.beta
        lam = args.lam
        loss = args.loss
        encoder_architecture = args.encoder_architecture
        decoder_architecture =args.decoder_architecture
        encoder_batch_norm = args.encoder_batch_norm
        verbose = args.verbose
        output_dir = args.output_dir
        return data_file,train_file,test_file,latent_dim, epochs, batch_size,optimizer,learning_rate,epsilon_std,beta,lam,loss,encoder_architecture,decoder_architecture,encoder_batch_norm,verbose,output_dir

    # read in the training and testing files
    def read_input_files(data_file, train_file,test_file):
        data = pd.read_csv(data_file, index_col=0)
        train = pd.read_csv(train_file, index_col = 0)
        test = pd.read_csv(test_file, index_col = 0)
        return data,train,test

    # find the training shape
    def get_input_dim(train):
        input_dim=train.shape[1]
        return input_dim



    # output training plot
    def train_plot(history,output_dir):
        plt.figure(figsize=(10, 5))
        plt.plot(history["loss"], label="Training data")
        plt.plot(history["val_loss"], label="Validation data")
        plt.ylabel("MSE + KL Divergence")
        plt.xlabel("No. Epoch")
        plt.title("VAE")
        plt.legend()
        #plt.show()
        plt.savefig(output_dir +'training_plot.png')



    # wrapper function to train autoencoder
    parser = define_arguments()

    data_file,train_file,test_file,latent_dim, epochs, batch_size,optimizer,learning_rate,epsilon_std,beta,lam,loss,encoder_architecture,decoder_architecture,encoder_batch_norm,verbose,output_dir = generate_arguments(parser)

    data, train,test = read_input_files(data_file,train_file, test_file)
  
    input_dim = get_input_dim(train)
  

    vae_model = VAE(
        input_dim=input_dim,
        latent_dim=latent_dim,
        batch_size=batch_size,
        encoder_batch_norm=encoder_batch_norm,
        epochs=epochs,
        learning_rate=learning_rate,
        encoder_architecture=encoder_architecture,
        decoder_architecture=decoder_architecture,
        beta=beta,
        lam=lam,
        verbose=verbose,
    )

    vae_model.compile_vae()
    vae_model.train(x_train=train,x_test=test)
    history = pd.DataFrame(vae_model.vae.history.history)
    train_plot(history, output_dir)
    vae_model.vae.evaluate(test)

    decoder= vae_model.decoder_block["decoder"]
    encoder= vae_model.encoder_block["encoder"]

    decoder.save_weights(output_dir, "decoder_weights.h5")
    encoder.save_weights(output_dir, "encoder_weights.h5")

    latent = pd.DataFrame(np.array(encoder.predict(data)[2]))
    latent.to_csv(output_dir+ "latent_data.csv")

    #train_plot(history, output_dir)
    
if __name__ == '__main__':
    main()
