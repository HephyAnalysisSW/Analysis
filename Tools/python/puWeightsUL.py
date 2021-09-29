''' PU reweighting from central POG json for UL v2
'''
#Standard imports
import ROOT, os
from correctionlib import _core

# Logging
import logging
logger = logging.getLogger(__name__)
    
puJson ="$CMSSW_BASE/src/Analysis/Tools/data/puReweightingData/LUM/"

def getPUReweight (year ="2016postVFP_UL"):
	fileName = puJson + year + "puWeights.json"
	logger.info( "Loaded 'pileup' from central POG json file %s", fileName )
	evaluator = _core.CorrectionSet.from_file(fileName)


