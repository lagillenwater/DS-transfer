from utils import *
import pandas as pd

parser = define_arguments()
data_file1,train_file1,test_file1, data_file2,train_file2,test_file2,latent_dim, epochs, batch_size,output_dir = generate_arguments(parser)
vae1,encoder1,decoder1,history1 = vanilla_vae(train_file1, test_file1, latent_dim,epochs, batch_size)
train_plot(history1,train_file1,output_dir)
vae2,encoder2,decoder2,history2 = vanilla_vae(train_file2, test_file2, latent_dim,epochs, batch_size)
train_plot(history2,train_file2,output_dir)

data1 = pd.read_csv(data_file1, index_col = 0)
data2 = pd.read_csv(data_file2, index_col = 0)

latent1 = encoder1.predict(data1, batch_size = batch_size)[2]
latent2 = encoder2.predict(data2, batch_size = batch_size)[2]

predict11 = pd.DataFrame(decoder1.predict(latent1))
predict11.columns = data1.columns
predict11.index = data1.index

predict22 = pd.DataFrame(decoder2.predict(latent2))
predict22.columns = data1.columns
predict22.index = data2.index

cross12 = pd.DataFrame(decoder1.predict(latent2))
cross12.columns = data2.columns
cross12.index = data2.index

cross21 = pd.DataFrame(decoder2.predict(latent1))
cross21.columns = data1.columns
cross21.index = data1.index

predict11.to_csv(output_dir+"predict11.csv")
predict22.to_csv(output_dir+"predict22.csv")
cross12.to_csv(output_dir+"cross12.csv")
cross21.to_csv(output_dir+"cross21.csv")
