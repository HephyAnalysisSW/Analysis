''' PU reweighting from central POG json for UL v2
'''
import os
import ROOT
from correctionlib import _core

# Logging
import logging
logger = logging.getLogger(__name__)
    
#puJson = "/users/priya.hussain/private/CMSSW_10_6_25/src/Analysis/Tools/data/puReweightingData/LUM/2016preVFP_UL/puWeights.json"
puJson = "/users/priya.hussain/private/CMSSW_10_6_25/src/Analysis/Tools/data/puReweightingData/LUM/"


def getPUReweight (nTrueInt, year ="UL2016",  weight="nominal" ):
	if "2016_preVFP" in year:
		tag = "2016preVFP_UL"
	elif "2016" in year:
		tag = "2016postVFP_UL"
	elif "2017" in year:
		tag = "2017_UL"
	else:
		tag = "2018_UL"
	fileName = puJson + tag + os.sep+ "puWeights.json"

	logger.debug( "Loaded 'pileup' from central POG json file %s", fileName )
	evaluator = _core.CorrectionSet.from_file(fileName)

	if "2016" in year:
		correction = evaluator["Collisions16_UltraLegacy_goldenJSON"]

	elif "2017" in year:
		correction = evaluator["Collisions17_UltraLegacy_goldenJSON"]

	else:
		correction = evaluator["Collisions18_UltraLegacy_goldenJSON"]

	PU = correction.evaluate(nTrueInt, weight) 

	return PU

