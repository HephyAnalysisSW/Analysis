#!/usr/bin/env python

# General
import os, sys
import ROOT

# Analysis
from Analysis.Tools.helpers import getVarValue, getObjDict

# RootTools
from RootTools.core.standard import *


import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel', action='store', nargs='?',  choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],   default='INFO', help="Log level for logging" )
argParser.add_argument('--sample',             action='store', type=str)
argParser.add_argument('--config_module',         action='store', type=str, default = "Analysis.MVA.cfg_examples", help = "config directory")
argParser.add_argument('--config',             action='store', type=str, default = "ttZ_dy_example", help="config")
argParser.add_argument('--output_directory',   action='store', type=str, default='.')
argParser.add_argument('--small',              action='store_true')

args = argParser.parse_args()

#Logger
import Analysis.Tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )

subDir = args.config

# MVA configuration
import importlib
configs = importlib.import_module(args.config_module)
config  = getattr( configs, args.config)

sample_names = []
found = False
for sample in config.training_samples:
    if args.sample == sample.name:
        found = True
        break # found it
    else:
        sample_names.append( sample.name )

if not found:
    logger.error( "Need sample to be one of %s, got %s", ",".join( sample_names ), args.sample )
    sys.exit()     

logger.info( "Processing sample %s", sample.name )

count  = int(sample.getYieldFromDraw( weightString="(1)" )["val"])
logger.info( "Found %i events for sample %s", count, sample.name )

if args.small:
    sample.reduceFiles(to=5)
    subDir += '_small'

# selection
if hasattr( config, "selectionString"):
    sample.addSelectionString( config.selectionString )
    logger.info( "Add selectionstring %s", config.selectionString )
else:
    logger.info( "Do not use selectionstring" )

if hasattr(sample, "selectionString"):
    logger.info( "Sample has %s", sample.selectionString) 

# where the output goes
output_file  = os.path.join( args.output_directory, "MVA-training", subDir, sample.name, sample.name + ".root" )

# reader
reader = sample.treeReader( \
    #variables = map( TreeVariable.fromString, config.read_variables),
    variables = config.read_variables + ( sample.read_variables if hasattr( sample, "read_variables") else []),
    sequence  = config.sequence,
    )

def fill_vector_collection( event, collection_name, collection_varnames, objects, nMax = 100):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects[:nMax]):
        for var in collection_varnames:
            if var in list(obj.keys()):
                if type(obj[var]) == type("string"):
                    obj[var] = int(ord(obj[var]))
                if type(obj[var]) == type(True):
                    obj[var] = int(obj[var])
                getattr(event, collection_name+"_"+var)[i_obj] = obj[var]

#filler
def filler( event ):

    r = reader.event

    # copy scalar variables
    for name, func in config.all_mva_variables.items():
        setattr( event, name, func(r, sample) )

    # copy vector variables
    for name, vector_var in config.mva_vector_variables.items():
        objs = vector_var["func"]( r, sample) 
        #print name, objs, vector_var['varnames']
        fill_vector_collection( event, name, vector_var['varnames'], objs, nMax = vector_var['nMax'] if 'nMax' in vector_var else None )

# Create a maker. Maker class will be compiled. 

# scalar variables
mva_variables = ["%s/F"%var for var in list(config.all_mva_variables.keys())]

# vector variables, if any
for name, vector_var in config.mva_vector_variables.items():
    #mva_variables.append( VectorTreeVariable.fromString(name+'['+','.join(vector_var['vars'])+']') )
    mva_variables.append ( VectorTreeVariable.fromString(name+'['+','.join(vector_var['vars'])+']', nMax = vector_var['nMax'] if 'nMax' in vector_var else None ) )
## FIs
#if hasattr( config, "FIs"):
#    FI_variables = ["FI_%s/F"%var for var in config.FIs.keys() ]
#else:
#    FI_variables = [] 

tmp_dir     = ROOT.gDirectory

dirname = os.path.dirname(output_file)
if not os.path.exists(dirname):
    os.makedirs(dirname)

outputfile = ROOT.TFile.Open(output_file, 'recreate')

# keep branches
if hasattr( config, "keep_branches" ):
    clonedTree = reader.cloneTree( config.keep_branches, newTreename = "Events", rootfile = outputfile )
else:
    clonedTree = None

outputfile.cd()
maker = TreeMaker(
    sequence  = [ filler ],
    variables = [ TreeVariable.fromString(var) if type(var)==type("") else var for var in  mva_variables], 
    treeName  = "Events",
    )

if clonedTree is not None:
    maker = maker.cloneWithoutCompile( externalTree = clonedTree )

tmp_dir.cd()

reader.start()
maker.start()

logger.info( "Starting event loop" )
counter=0
while reader.run():

    maker.run()
    counter += 1
    if counter%10000 == 0:
        logger.info("Written %i events.", counter)

nEventsTotal = maker.tree.GetEntries()

maker.tree.Write()
outputfile.Close()
logger.info( "Written %s", output_file)
#
#      # Destroy the TTree
maker.clear()
logger.info( "Written %i events to %s",  nEventsTotal, output_file )
