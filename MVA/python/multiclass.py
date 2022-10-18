#!/usr/bin/env python

import ROOT, os
from RootTools.core.standard import *
import Analysis.Tools.syncer as syncer
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('/tmp/delete.png')

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--config_module',      action='store', type=str, default = "Analysis.MVA.cfg_examples", help = "config directory")
argParser.add_argument('--config',             action='store', type=str, default = "ttZ_dy_example", help="config")
argParser.add_argument('--name',               action='store', type=str,   default='default', help="Name of the training")
argParser.add_argument('--variable_set',       action='store', type=str,   default='mva_variables', help="List of variables for training")
argParser.add_argument('--output_directory',   action='store', type=str,   default=os.path.expandvars('/groups/hephy/cms/$USER/'))
argParser.add_argument('--input_directory',    action='store', type=str,   default=os.path.expandvars("/scratch-cbe/users/$USER/example-ntuples-ttZ_dy/MVA-training"))
argParser.add_argument('--small',              action='store_true', help="small?")
argParser.add_argument('--add_LSTM',           action='store_true', help="add LSTM?")

args = argParser.parse_args()

if args.add_LSTM:   args.name+="_LSTM"
input_postfix = ''
if args.small:
    args.name    +="_small"
    input_postfix = '_small'

#Logger
import Analysis.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )

# MVA configuration
import importlib
configs = importlib.import_module(args.config_module)
config  = getattr( configs, args.config)

import uproot
import awkward
import numpy as np
import pandas as pd
#import h5py

#########################################################################################
# variable definitions

import Analysis.Tools.user as user

# directories
plot_directory   = os.path.join( user. plot_directory, 'MVA', args.name, args.config )
output_directory = os.path.join( args.output_directory, 'models', args.name, args.config)

# fix random seed for reproducibility
np.random.seed(1)

# get the training variable names
mva_variables = [ mva_variable[0] for mva_variable in getattr(config, args.variable_set) ]

n_var_flat   = len(mva_variables)

df_file = {}
for i_training_sample, training_sample in enumerate(config.training_samples):
    upfile_name = os.path.join(os.path.expandvars(args.input_directory), args.config+input_postfix, training_sample.name, training_sample.name+'.root')
    logger.info( "Loading upfile %i: %s from %s", i_training_sample, training_sample.name, upfile_name)
    upfile = uproot.open(upfile_name)
    df_file[training_sample.name]  = upfile["Events"].pandas.df(branches = mva_variables )
    # enumerate
    df_file[training_sample.name]['signal_type'] =  np.ones(len(df_file[training_sample.name])) * i_training_sample

df = pd.concat([df_file[training_sample.name] for training_sample in config.training_samples])

#df = df.dropna() # removes all Events with nan -> amounts to M3 cut

# split dataset into Input and output data
dataset = df.values

# number of samples with 'small'
n_small_samples = 10000

# small
if args.small:
    dataset = dataset[:n_small_samples]

X  = dataset[:,0:n_var_flat]

# regress FI
Y = dataset[:, n_var_flat]

from sklearn.preprocessing import label_binarize
classes = list(range(len(config.training_samples)))

if len(config.training_samples) == 2: Y = label_binarize(Y, classes=classes+[-1])[:,:2]
else: Y = label_binarize(Y, classes=classes)

# loading vector branches for LSTM
if args.add_LSTM:
    vector_branches = ["mva_Jet_%s" % varname for varname in config.lstm_jetVarNames]
    max_timestep = config.lstm_jets_maxN # for LSTM

    vec_br_f  = {}

    for i_training_sample, training_sample in enumerate(config.training_samples):
        upfile_name = os.path.join(os.path.expandvars(args.input_directory), args.config+input_postfix, training_sample.name, training_sample.name+'.root')
        logger.info( "Loading vector branches %i: %s from %s", i_training_sample, training_sample.name, upfile_name)
        with uproot.open(upfile_name) as upfile:
            vec_br_f[i_training_sample]   = {}
            for name, branch in upfile["Events"].arrays(vector_branches).items():
                vec_br_f[i_training_sample][name] = branch.pad(max_timestep)[:,:max_timestep].fillna(0)

    vec_br = {name: awkward.JaggedArray.concatenate( [vec_br_f[i_training_sample][name] for i_training_sample in range(len(config.training_samples))] ) for name in vector_branches}
    if args.small:
        for key, branch in vec_br.items():
            vec_br[key] = branch[:n_small_samples]

    # put columns side by side and transpose the innermost two axis
    len_samples = len(list(vec_br.values())[0])
    V           = np.column_stack( [vec_br[name] for name in vector_branches] ).reshape( len_samples, len(vector_branches), max_timestep).transpose((0,2,1))

# split data into train and test, test_size = 0.2 is quite standard for this
from sklearn.model_selection import train_test_split

options = {'test_size':0.2, 'random_state':7, 'shuffle':True}
if args.add_LSTM:
    X_train, X_test, Y_train, Y_test, V_train, V_test = train_test_split(X, Y, V, **options)
    validation_data = ( [X_test,  V_test], Y_test )
    training_data   =   [X_train, V_train]
else:
    X_train, X_test, Y_train, Y_test                  = train_test_split(X, Y, **options)
    validation_data = ( X_test,  Y_test)
    training_data   =   X_train

#########################################################################################
# define model (neural network)
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten, LSTM, Concatenate
from keras.layers import BatchNormalization
from keras.utils import np_utils

# flat layers
flat_inputs = Input(shape=(n_var_flat, ))
x = BatchNormalization(input_shape=(n_var_flat, ))(flat_inputs)
x = Dense(n_var_flat*2, activation='sigmoid')(x)
x = Dense(n_var_flat+5, activation='sigmoid')(x)

inputs = flat_inputs

# LSTMs
if args.add_LSTM:
    vec_inputs = Input(shape=(max_timestep, len(vector_branches),) )
    v = LSTM(10, activation='relu', input_shape=(max_timestep, len(vector_branches)))( vec_inputs )
    x = Concatenate()( [x, v])
    inputs = ( flat_inputs, vec_inputs)

outputs = Dense(len(config.training_samples), kernel_initializer='normal', activation='sigmoid')(x)
model = Model( inputs, outputs )

#model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mean_absolute_percentage_error'])
model.summary()

# define callback for early stopping
import tensorflow as tf
callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=3) # patience can be higher if a more accurate result is preferred
                                                                        # I would recommmend at least 3, otherwise it might cancel too early

# train the model
batch_size = 1024*6
history = model.fit(training_data,
                    Y_train,
                    sample_weight = None,
                    epochs=500,
                    batch_size=batch_size,
                    #verbose=0, # switch to 1 for more verbosity, 'silences' the output
                    #validation_split=0.1
                    validation_data = validation_data,
                    callbacks=[callback],
                   )
print('training finished')

# plot the metrics
from plot_metrics import plot_metrics
plot_metrics(history, plot_directory)

# saving
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

output_file = os.path.join(output_directory, 'multiclass_model.h5')
model.save(output_file)
logger.info("Written model to: %s", output_file)
