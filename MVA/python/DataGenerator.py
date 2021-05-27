#!/usr/bin/env python

import os
from RootTools.core.standard import *

# Logging
import logging
logger = logging.getLogger(__name__)

import uproot
import awkward
import numpy as np
import pandas as pd

def chunk( tot, n_split, index):
    ''' Implements split of number into n_split chunks
        https://math.stackexchange.com/questions/2975936/split-a-number-into-n-numbers
        Return tuple (start, stop)
    '''
        
    d = int( tot / n_split )
    r = tot % n_split

    #return [d+1]*r + [d]*(n_split-r)
    if index<r:
        return ( (d+1)*index, (d+1)*(index+1) )
    else:
        return ( (d+1)*r + d*(index-r), (d+1)*r + d*(index-r+1) )


from keras.utils import Sequence
class DataGenerator(Sequence):

    default_options = {'test_size':0.2, 'random_state':42, 'shuffle':True}

    def __init__( self, 
            config, 
            n_split, 
            input_directory,
            return_test_data=False, 
            variable_set = "mva_variables",
            jet_LSTM     = False,
            small        = False,
            **kwargs
                ):
        '''
        DataGenerator for the keras training framework.
        config:     A configuration class.
        n_split:    Number of chunks the data should be split into.
        input_directory:    Location of the training data.
        return_test_data:   If True, the test data is returned. Otherwise, the training data is returned.
        variable_set:   Name of the attribute of the confige that has the variables.
        jet_LSTM:   Add jet data for an LSTM?
        small:      Process only small subset (10k events)?
        kwargs:     Additional configuration for the train_test_split method
        '''

        # Into how many chunks we split
        self.n_split        = n_split
        if not n_split>0 or not type(n_split)==int:
            raise RuntimeError( "Need to split in positive integer. Got %r"%n )

        # input_directory, i.e., the output of make_ntuple.py
        self.input_directory= input_directory

        # MVA config
        self.config  = config

        # train_test_split options
        self.options = DataGenerator.default_options
        self.options.update(kwargs)
    
        # variables to return 
        self.variable_set   = variable_set
        self.mva_variables  = [ mva_variable[0] for mva_variable in getattr(self.config, variable_set) ]

        # recall the index we loaded
        self.index          = -1

        # whether or not to return the test data
        self.return_test_data = return_test_data

        # options
        self.small          = small
        self.jet_LSTM       = jet_LSTM

        # load the files (don't load the data)
        self.load_files()
    
    def validation_data_generator( self ):
        ''' Make a new DataGenerator that returns the test data'''
        return DataGenerator( self.config, self.n_split, self.input_directory, 
                return_test_data = True,
                variable_set = self.variable_set,
                jet_LSTM     = self.jet_LSTM,
                small        = self.small,
                **self.options
        )

    # interface to Keras
    def __len__( self ):
        return self.n_split

    def load_files( self ):
        logger.info( "Initializing " + ("test data." if self.return_test_data else "training data."))   

        self.n_total     = {}
        self.events_tree = {}
        for i_training_sample, training_sample in enumerate(self.config.training_samples):
            upfile_name = os.path.join(os.path.expandvars(self.input_directory), self.config.__name__+('_small' if self.small else ''), training_sample.name, training_sample.name+'.root')
            upfile_     = uproot.open( upfile_name )
            logger.info( "upfile %i: %s from %s", i_training_sample, training_sample.name, upfile_name)

            # total number of events in the tree
            self.n_total[training_sample.name]      = len(upfile_["Events"])
            self.events_tree[training_sample.name]  = upfile_["Events"]

    def load( self, index ):
        logger.info( "Loading data batch %i/%i (%s)" % (index, self.n_split, ("test" if self.return_test_data else "train")) )    

        df_file = {}
        for i_training_sample, training_sample in enumerate(self.config.training_samples):

            # total number of events in the tree
            entrystart, entrystop = chunk( self.n_total[training_sample.name], self.n_split, index ) 

            df_file[training_sample.name]  = self.events_tree[training_sample.name].pandas.df(branches = self.mva_variables, entrystart=entrystart, entrystop=entrystop)
            # enumerate
            df_file[training_sample.name]['signal_type'] =  np.ones(len(df_file[training_sample.name])) * i_training_sample

        df = pd.concat([df_file[training_sample.name] for training_sample in self.config.training_samples])

        # split dataset into Input and output data
        dataset = df.values

        # number of samples with 'small'
        n_small_samples = 10000

        # small
        if self.small:
            dataset = dataset[:n_small_samples]

        n_var_flat     = len(self.mva_variables)
        X  = dataset[:,0:n_var_flat]

        Y = dataset[:, n_var_flat] 

        from sklearn.preprocessing import label_binarize
        classes = range(len(self.config.training_samples))

        if len(self.config.training_samples) == 2: Y = label_binarize(Y, classes=classes+[-1])[:,:2]
        else: Y = label_binarize(Y, classes=classes)

        # loading vector branches for LSTM
        if self.jet_LSTM:
            vector_branches = ["mva_Jet_%s" % varname for varname in self.config.lstm_jetVarNames]
            max_timestep = self.config.lstm_jets_maxN # for LSTM

            vec_br_f = {}

            for i_training_sample, training_sample in enumerate(self.config.training_samples):
                logger.info( "Loading vector branches %i from %s (%s)", i_training_sample, training_sample.name, ("test data." if self.return_test_data else "training data."))
                vec_br_f[i_training_sample] = {}

                entrystart, entrystop = chunk( self.n_total[training_sample.name], self.n_split, index ) 

                for name, branch in self.events_tree[training_sample.name].arrays(vector_branches, entrystart=entrystart, entrystop=entrystop).iteritems():
                    vec_br_f[i_training_sample][name] = branch.pad(max_timestep)[:,:max_timestep].fillna(0)

            vec_br = {name: awkward.JaggedArray.concatenate( [vec_br_f[i_training_sample][name] for i_training_sample in range(len(self.config.training_samples))] ) for name in vector_branches}
            if self.small:
                for key, branch in vec_br.iteritems():
                    vec_br[key] = branch[:n_small_samples]

            # put columns side by side and transpose the innermost two axis
            len_samples = len(vec_br.values()[0])
            V           = np.column_stack( [vec_br[name] for name in vector_branches] ).reshape( len_samples, len(vector_branches), max_timestep).transpose((0,2,1))

        # split data into train and test, test_size = 0.2 is quite standard for this
        from sklearn.model_selection import train_test_split

        if self.jet_LSTM:
            X_train, X_test, Y_train, Y_test, V_train, V_test = train_test_split(X, Y, V, **self.options)
            if self.return_test_data:
                self.data = ( [X_test,  V_test],  Y_test )
            else:
                self.data = ( [X_train, V_train], Y_train )
        else:
            X_train, X_test, Y_train, Y_test                  = train_test_split(X, Y, **self.options)
            if self.return_test_data:
                self.data = ( X_test,  Y_test)
            else:
                self.data   = ( X_train, Y_train)

        # recall what index we just loaded
        self.index = index

    def __getitem__(self, index):
        if index == self.index:
            return self.data
        else:
            self.load( index )
            return self.data 

if __name__=='__main__':

    import TMB.MVA.configs as configs
    data = DataGenerator(
        config  = configs.ttZ_dy,
        n_split = 1000,
        jet_LSTM = True,
        input_directory = os.path.expandvars("/scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/"),
        ) 
