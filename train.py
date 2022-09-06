from vae import VAE
import argparse
import pandas as pd

#Define arguments for each required and optional input
def define_arguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ## Required inputs
    parser.add_argument("--train-file",dest="train_file",required=True,help="training dataset")
    parser.add_argument("--test-file",dest="test_file",required=True,help="testing dataset")
    parser.add_argument("--encoder-architecture",dest="encoder_architecture",required=True,help="encoder architecture")
    parser.add_argument("--decoder-architecture",dest="decoder_architecture",required=True,help="decoder architecture")

    # optional inputs
    parser.add_argument("--epochs", dest = "epochs", require=False, help="epochs", default=50)
    parser.add_argument("--latent-dim",dest="latent_dim",required=False,help="latent dimensions",default = 16)
    parser.add_argument("--batch-size",dest="bat_size",required=False,default=128,help="batch size")
    parser.add_argument("--optimizer",dest="optimizer",required=False,default='adam',help="optimizer")
    parser.add_argument("--learning-rate",dest="learning_rate",required=False,default=0.0005,help="learning_rate")
    parser.add_argument("--epsilon",dest="epsilon_std",required=False,default=1.0,help="epsilon")
    parser.add_argument("--lambda",dest="lam",required=False,default=0.0,help="lambda")
    parser.add_argument("--beta",dest="beta",required=False,default=1.0,help="beta")
    parser.add_argument("--loss",dest="loss",required=False,default="binary_crossentropy",help="loss function")
    parser.add_argument("--encoder_batch_norm",dest="encoder_batch_norm",required=False,default=True,help="encoder batch norm")
    parser.add_argument("--verbose",dest="verbose",required=False,default=True,help="verbose")
    return parser

# generate arguments
def generate_arguments(parser):

    #Generate argument parser and define arguments
    parser = define_arguments()
    args = parser.parse_args()

    train_file = args.train_file
    test_file = args.test_file
    latent_dim = args.latent_dim
    epochs = args.epochs
    batch_size = args.batch_size
    optimizer = args.optimizer
    learning_rate = args.learning_rate
    epsilon_std = args.epsilon_std
    beta =K.variable(args.beta)
    lam = K.variable(args.lam)
    loss = args.loss
    encoder_architecture = args.encoder_architecture
    decoder_architecture = args.decoder_architecture
    encoder_batch_norm = args.encoder_batch_norm
    verbose = args.verbose
    return train,test,latent_dim, epochs, batch_size,optimizer,learning_rate,epsilon_std,beta,lam,loss,encoder_architecture,decoder_architecture,encoder_batch_norm,verbose

def read_input_files(train_file,test_file):
    train = pd.read_csv(train_file)
    test = pd.read_csv(test_file)
    return train,test
