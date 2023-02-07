import pandas as pd
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
import argparse

def main():

    def define_arguments():
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        # Required input
        parser.add_argument("--i", dest = "input_file", required = True, help = "input csv file")
        parser.add_argument("--o", dest = "output_dir", required = True, help = "output_dir")
        return parser

    def generate_arguments(parser):
        # Generate argument parser and define arguments
        parser=define_arguments()
        args = parser.parse_args()

        input_file = args.input_file
        output_dir = args.output_dir
        return input_file, output_dir

    # load inputs
    def load_input(input_file):
        # reading in the data
        data = pd.read_csv(input_file, index_col = 0)
        return data

    # scale data
    def min_max_scale(data):
        # scale data prior to input into the model
        scaler = preprocessing.MinMaxScaler(feature_range=(0,1))
        scaler.fit(data)
        minmax_data = pd.DataFrame(scaler.transform(data))
        minmax_data.columns= data.columns
        minmax_data.index = data.index
        return(minmax_data)

    def split_data(scaled_data):
        
        # split to training and test data
        # parameterize in the future
        test_split = .3
        seed = 42

        # Split dataxb
        train_df, test_df = train_test_split(
        scaled_data,
        test_size=test_split,
        random_state=seed,
        #    stratify = sample_annotations.cell_type,
        )
        return train_df,test_df

    def output_train_test(train_df, test_df, output_dir, input_file):
        train_output = output_dir+ "train_" + input_file.split('/')[-1]
        train_df.to_csv(train_output)
        test_output = output_dir+ "test_" + input_file.split('/')[-1]
        test_df.to_csv(test_output)

        
    def process_wrapper():
        parser = define_arguments()
        input_file, output_dir = generate_arguments(parser)
        data = load_input(input_file)
        scaled_data = min_max_scale(data)
        train_df, test_df = split_data(scaled_data)
        output_train_test(train_df, test_df, output_dir, input_file)

    process_wrapper()

if __name__ == '__main__':
    main()



