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
argParser.add_argument('--split',              action='store', type=int,   default=1, help="In how many batches do you want tos split the training data?")
argParser.add_argument('--variable_set',       action='store', type=str,   default='mva_variables', help="List of variables for training")
argParser.add_argument('--output_directory',   action='store', type=str,   default=os.path.expandvars('/mnt/hephy/cms/$USER/models/'))
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
# Training data 

import Analysis.Tools.user as user 

# directories
plot_directory   = os.path.join( user. plot_directory, 'MVA', args.name, args.config )
output_directory = os.path.join( args.output_directory, args.name, args.config) 

# fix random seed for reproducibility
np.random.seed(1)

from Analysis.MVA.DataGenerator import DataGenerator 
training_data_generator = DataGenerator(
        config  = config,
        n_split = 1000,
        jet_LSTM = args.add_LSTM,
        input_directory = os.path.expandvars(args.input_directory),
        small = args.small
        )
validation_data_generator = training_data_generator.validation_data_generator()

#########################################################################################
# define model (neural network)
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten, LSTM, Concatenate
from keras.layers import BatchNormalization
from keras.utils import np_utils


# flat layers

n_var_flat = len(getattr(config, args.variable_set))
flat_inputs = Input(shape=(n_var_flat, ))
x = BatchNormalization(input_shape=(n_var_flat, ))(flat_inputs)
x = Dense(n_var_flat*2, activation='sigmoid')(x)
x = Dense(n_var_flat+5, activation='sigmoid')(x)

inputs = flat_inputs

# LSTMs
if args.add_LSTM:
    max_timestep = config.lstm_jets_maxN
    vector_branches = ["mva_Jet_%s" % varname for varname in config.lstm_jetVarNames]

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
history = model.fit_generator(
                    training_data_generator,
                    epochs=500, 
                    #verbose=0, # switch to 1 for more verbosity, 'silences' the output
                    #validation_split=0.1
                    validation_data = validation_data_generator,
                    callbacks=[callback],
                    shuffle=False
                   )
print('training finished')

# saving
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

output_file = os.path.join(output_directory, 'multiclass_model.h5')
model.save(output_file)
logger.info("Written model to: %s", output_file)
