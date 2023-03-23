'''Implementation of b-tagging reshaping
'''

#Standard Imports
import ROOT
import os, sys, subprocess, shutil, uuid
import numpy as np
import awkward as ak
from correctionlib import _core
#Logging
import logging
logger = logging.getLogger(__name__)

#RootTools Imports
from RootTools.core.standard import *

def toFlavourKey(pdgId):
    if abs(pdgId)==5: return 0
    if abs(pdgId)==4: return 1
    return 2

flavourSys = {
    5:{'central', 'up_jes', 'down_jes', 'up_lf', 'down_lf', 'up_hf', 'down_hf', 'up_hfstats1', 'down_hfstats1', 'up_hfstats2', 'down_hfstats2', 'up_lfstats1', 'down_lfstats1', 'up_lfstats2', 'down_lfstats2'},
    4:{'central', 'up_cferr1', 'down_cferr1', 'up_cferr2', 'down_cferr2'},
    0:{'central', 'up_jes', 'down_jes', 'up_lf', 'down_lf', 'up_hf', 'down_hf', 'up_hfstats1', 'down_hfstats1', 'up_hfstats2', 'down_hfstats2', 'up_lfstats1', 'down_lfstats1', 'up_lfstats2', 'down_lfstats2'},
}



class BTagReshaping:

    def __init__(self, fastSim=False, year='UL2016', tagger='DeepJet' ):

        if year not in ['UL2016', 'UL2016_preVFP', 'UL2017', 'UL2018']:
            raise Exception("b tag SF for year %s not known"%year)

        self.fastSim = fastSim  #Whether or not FS SF are to be used
        self.year = year
        self.tagger = tagger
        self.dataDir = "$CMSSW_BASE/src/Analysis/Tools/data/btagEfficiencyData/BTV"

        #Input files
        b_file = "btagging.json"
        #c_file = "ctagging.json"

        if year == 'UL2016_preVFP':
           tag = '2016preVFP_UL'
           self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, b_file ) )
    	if year == 'UL2016':
    	   tag = '2016postVFP_UL'
           self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, b_file ) )
        if year == 'UL2017':
    	   tag = '2017_UL'
           self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, b_file ) )
        if year == 'UL2018':
    	   tag = '2018_UL'
           self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, b_file ) )


        logger.info ( "Loading scale factors from %s", self.scaleFactorFile )
    	self.correction = _core.CorrectionSet.from_file(self.scaleFactorFile)


    def getbtagSF(self, j):
        # #Create a dictionary that fetches the right values
        j['jetSF'] = {var: self.correction["deepJet_shape"].evaluate(var, j['hadronFlavour'], abs(j['eta']), min(10000.,j['pt']), j['btagDeepFlavB']) for var in list(flavourSys[abs(j['hadronFlavour'])])}

        
