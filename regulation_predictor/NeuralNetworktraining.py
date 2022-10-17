import tensorflow as tf
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split
from collections import Counter
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, InputLayer, BatchNormalization, Dropout
from sklearn.metrics import classification_report
from tensorflow.keras.utils import to_categorical
from sklearn.compose import ColumnTransformer
from tensorflow.keras import regularizers
import os

workdir = os.getcwd()

NNdata = pd.read_csv(f'{workdir}/Data_Sets/trainingsets.csv', index_col=0)
NNdata.replace([np.inf, -np.inf], np.nan, inplace=True)
NNdata.dropna(inplace=True)

y = NNdata[['is_upreg']]
x = NNdata.drop('is_upreg',axis=1)

for k in [2]:
    for col in x.columns:
        x[col+'_'+str(k)] = x[col].apply(lambda x: x**k)


x_train, x_test, y_train, y_test = train_test_split(x,y,random_state= 118, test_size= .17, stratify=y)
ct = ColumnTransformer([('numeric',StandardScaler(),x.columns)])
x_train = ct.fit_transform(x_train)
x_test = ct.fit_transform(x_test)
le = LabelEncoder()
y_train = le.fit_transform(y_train.astype('str'))
y_test = le.fit_transform(y_test.astype('str'))
y_train = to_categorical(y_train)
y_test = to_categorical(y_test)


#Hyperparameter tunning
learning_rate = 0.01
batchsize = 22
nepoc = 300

model = Sequential(name = 'sRNApredicter')
model.add(InputLayer(input_shape=(x_train.shape[1])))
model.add(BatchNormalization(axis=-1,momentum=0.99,epsilon=0.001,center=True,scale=True,beta_initializer="zeros",gamma_initializer="ones",
    moving_mean_initializer="zeros",moving_variance_initializer="ones",))
model.add(Dense(128,kernel_regularizer=regularizers.L1L2(l1=1e-5, l2=1e-4),bias_regularizer=regularizers.L1(1e-3),
    activity_regularizer=regularizers.L2(1e-5), activation = 'relu'))
model.add(Dense(128,kernel_regularizer=regularizers.L1L2(l1=1e-5, l2=1e-4),bias_regularizer=regularizers.L1(1e-4),
    activity_regularizer=regularizers.L2(1e-5), activation = 'relu'))
model.add(Dense(128,kernel_regularizer=regularizers.L1L2(l1=1e-5, l2=1e-4),bias_regularizer=regularizers.L1(1e-4),
    activity_regularizer=regularizers.L2(1e-5), activation = 'relu'))
model.add(Dropout(rate = .2, noise_shape=None, seed=None))
model.add(Dense(2,kernel_regularizer=regularizers.L1L2(l1=1e-5, l2=1e-4),bias_regularizer=regularizers.L1(1e-3),
    activity_regularizer=regularizers.L2(1e-5),activation = 'sigmoid'))

optimizer = tf.keras.optimizers.Adam(learning_rate = learning_rate)
model.compile(optimizer=optimizer,loss='categorical_crossentropy',metrics=['accuracy'])

model.fit(x_train, y_train, epochs= nepoc, batch_size = batchsize , verbose = 2, validation_split = 0.20)


loss, acc = model.evaluate(x_test, y_test, verbose=2)
print('sRNApredicter, accuracy: {:5.2f}%'.format(100 * acc))


model.save(f'{workdir}/regulation_predictor/sRNAprediction.h5', overwrite=True, include_optimizer=True,save_format=None, signatures=None, options=None, save_traces=True)

print('model saved as sRNAprediction')

