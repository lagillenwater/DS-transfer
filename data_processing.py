import pandas as pd
from sklearn import preprocessing
from sklearn.model_selection import train_test_split

# reading in the data
data = pd.read_csv("recount3_T21.csv", index_col = 0)

# scale data prior to input into the model
scaler = preprocessing.MinMaxScaler(feature_range=(0,1))
scaler.fit(data)
minmax_data = scaler.transform(data)
minmax_data.shape


# split to training and test data
test_split = .3
seed = 42

# Split dataxb
train_df, test_df = train_test_split(
    minmax_full_data,
    test_size=test_split,
    random_state=seed,
#    stratify = sample_annotations.cell_type,
)

test_df, valid_df = train_test_split(
    test_df,
    test_size=0.2,
    random_state=seed
)



print(train_df.shape)
print(test_df.shape)
print(valid_df.shape)





