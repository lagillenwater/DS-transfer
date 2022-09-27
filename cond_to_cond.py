import tensorflow as tf
from tensorflow import keras
import pandas as pd
from tensorflow.keras import backend as K

#def main():
#Define arguments for each required and optional input

def define_arguments():
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ## Required inputs
    parser.add_argument("--condition-1-directory",dest="cond1_dir",required=True,help="condition 1 directory")        
    parser.add_argument("--condition-2-directory",dest="cond2_dir",required=True,help="condition 2 directory")
    parser.add_argument("--output-dir",dest="output_dir",required=True,help="output directory")
    return parser

def generate_arguments(parser):
            #Generate argument parser and define arguments
    parser = define_arguments()
    args = parser.parse_args()

    cond1_dir = args.cond1_dir
    cond2_dir = args.cond2_dir
    output_dir = args.output_dir
    return cond1_dir,cond2_dir,output_dir

def lrelu(x, alpha=0.3):
    return tf.maximum(x, tf.multiply(x, alpha))

def load_models(cond1_dir,cond2_dir):
    cond1_model = keras.models.load_model(cond1_dir+"model.h5", custom_objects={'lrelu':lrelu, 'K':K})
    cond2_model = keras.models.load_model(cond2_dir+"model.h5", custom_objects={'lrelu':lrelu, 'K':K})
    cond1_decoder = cond1_model.load_weights(cond1_dir + "model_weights.h5")
    cond2_decoder = cond1_model.load_weights(cond2_dir + "model_weights.h5")
    return cond1_model,cond2_model, cond1_decoder, cond2_decoder

def load_latent_space(cond1_dir, cond2_dir):
    cond1_latent_space = pd.read_csv(cond1_dir+"latent_data.csv")
    cond2_latent_space = pd.read_csv(cond2_dir+"latent_data.csv")
    return cond1_latent_space, cond2_latent_space

def get_decoders(cond1_model, cond2_model):
    cond1_decoder= cond1_model.decoder_block["decoder"]
    cond2_decoder = cond2_model.decoder_block["decoder"]
    return cond1_decoder,cond2_decoder
