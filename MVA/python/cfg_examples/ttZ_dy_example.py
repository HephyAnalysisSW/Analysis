#!/usr/bin/env python

# Standard imports
from operator                   import attrgetter
from math                       import pi, sqrt, cosh, cos, acos
import ROOT, os

# RootTools
from RootTools.core.standard     import *

# helpers
from Analysis.Tools.helpers          import deltaPhi, deltaR2, deltaR, getCollection, getObjDict
from Analysis.Tools.WeightInfo       import WeightInfo

from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons

jetVars          = ['pt/F', 'eta/F', 'phi/F', 'btagDeepB/F', 'btagDeepFlavB/F', 'index/I']

jetVarNames      = [x.split('/')[0] for x in jetVars]

lstm_jets_maxN   = 10
lstm_jetVars     = ['pt/F', 'eta/F', 'phi/F', 'btagDeepFlavB/F', 'btagDeepFlavC/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'muEF/F', 'puId/F', 'qgl/F']
lstm_jetVarNames = [x.split('/')[0] for x in lstm_jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'mvaTOP/F', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
                    "nBTag/I",
                    "nJetGood/I",
                    "nlep/I",
                    "m3/F",
                    "JetGood[%s]"%(",".join(jetVars)),
                    "Jet[%s]"%(",".join(lstm_jetVars)),
                    "lep[%s]"%(",".join(lepVars)),
                    "met_pt/F", "met_phi/F",
                    "l1_pt/F",
                    "l1_eta/F",
                    "l1_phi/F",
                    "l2_pt/F",
                    "l2_eta/F",
                    "l2_phi/F",
                    #"l3_pt/F",
                    #"l3_eta/F",
                    #"l3_phi/F",
                    "l1_mvaTOP/F",
                    "l2_mvaTOP/F",
                    #"l3_mvaTOP/F",
                    "year/I",
                    ]
# sequence 
sequence = []

# Fisher informations
FIs = {
}

# Reco b-Jet Filter
def isBJet( j, tagger='DeepCSV', year=2016 ):
    if tagger == 'CSVv2':
        if year == 2016:
            # https://twiki.cern.ch/twikix/bin/viewauth/CMS/BtagRecommendation80XReReco
            return j['btagCSVV2'] > 0.8484 
        elif year == 2017:
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
            return j['btagCSVV2'] > 0.8838 
        elif year == 2018:
            # UPDATE WHEN AVAILABLE
            return j['btagCSVV2'] > 0.8838 
        else:
            raise NotImplementedError
    elif tagger == 'DeepCSV':
        if year == 2016:
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
            return j['btagDeepB'] > 0.6321
        elif year == 2017:
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
            return j['btagDeepB'] > 0.4941
        elif year == 2018:
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
            return j['btagDeepB'] > 0.4184
        else:
            raise NotImplementedError

def make_jets( event, sample ):
    event.jets     = [getObjDict(event, 'JetGood_', jetVarNames, i) for i in range(int(event.nJetGood))] 
    event.bJets    = [j for j in event.jets if isBJet(j, year=event.year) and abs(j['eta'])<=2.4]
sequence.append( make_jets )

def get_mll(event, sample):
    event.mll = sqrt(2*(event.l1_pt)*(event.l2_pt)*(cosh(event.l1_eta-event.l2_eta)-cos(event.l1_phi-event.l2_phi)))
sequence.append(get_mll)

def getDeltaR(event, sample):
    if len(event.jets) >= 1:
        event.minDRjet_l1       = min(deltaR({'eta':event.jets[j]['eta'],  'phi':event.jets[j]['phi']}, {'eta':event.l1_eta, 'phi':event.l1_phi})  for j in range(len(event.jets)) )
        event.minDRjet_l2       = min(deltaR({'eta':event.jets[j]['eta'],  'phi':event.jets[j]['phi']}, {'eta':event.l2_eta, 'phi':event.l2_phi})  for j in range(len(event.jets)))
    else: 
        event.minDRjet_l1      = -1 
        event.minDRjet_l2      = -1
    if len(event.bJets) >= 1:
        event.minDRbjet_l1      = min(deltaR({'eta':event.bJets[b]['eta'], 'phi':event.bJets[b]['phi']}, {'eta':event.l1_eta, 'phi':event.l1_phi}) for b in range(len(event.bJets) ))
        event.minDRbjet_l2      = min(deltaR({'eta':event.bJets[b]['eta'], 'phi':event.bJets[b]['phi']}, {'eta':event.l2_eta, 'phi':event.l2_phi}) for b in range(len(event.bJets )))
    else: 
        event.minDRbjet_l1      = -1 
        event.minDRbjet_l2      = -1
sequence.append(getDeltaR)

all_mva_variables = {

# global event properties     
     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
     "mva_nBTag"                 :(lambda event, sample: event.nBTag),
#     "mva_nlep"                  :(lambda event, sample: event.nlep),

     "mva_met_pt"                :(lambda event, sample: event.met_pt),
     "mva_l1_pt"                 :(lambda event, sample: event.l1_pt),
     "mva_l1_eta"                :(lambda event, sample: event.l1_eta),
     "mva_l2_pt"                 :(lambda event, sample: event.l2_pt),
     "mva_l2_eta"                :(lambda event, sample: event.l2_eta),
     "mva_l1_phi"                :(lambda event, sample: event.l1_phi),
     "mva_l2_phi"                :(lambda event, sample: event.l2_phi),

     "mva_mll"                   :(lambda event, sample: event.mll),
     "mva_l1_relIso"                :(lambda event, sample: event.lep_pfRelIso03_all[0]),
     "mva_l2_relIso"                :(lambda event, sample: event.lep_pfRelIso03_all[1]),

     "mva_ht"                    :(lambda event, sample: sum( [j['pt'] for j in event.jets] ) ),

     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0]   if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
     "mva_jet2_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=3 else 0),
     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -10),

     "mva_jet3_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=4 else 0),
     "mva_jet4_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=5 else 0),
     "mva_jet5_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=6 else 0),

     "mva_minDR_jl1"                :(lambda event, sample: event.minDRjet_l1  ),
     "mva_minDR_bjl1"               :(lambda event, sample: event.minDRbjet_l1 ),
     "mva_minDR_jl2"                :(lambda event, sample: event.minDRjet_l2  ),
     "mva_minDR_bjl1"               :(lambda event, sample: event.minDRbjet_l2 ),

                }

def lstm_jets(event, sample):
    jets = [ getObjDict( event, 'Jet_', lstm_jetVarNames, event.JetGood_index[i] ) for i in range(int(event.nJetGood)) ]
    return jets

# for the filler
mva_vector_variables    =   {
    "mva_Jet":  {"func":lstm_jets, "name":"Jet", "vars":lstm_jetVars, "varnames":lstm_jetVarNames}
}

## Using all variables
mva_variables_ = list(all_mva_variables.keys())
mva_variables_.sort()
mva_variables  = [ (key, value) for key, value in all_mva_variables.items() if key in mva_variables_ ]

# keep these branches in ntuple making
keep_branches = ["GenMET_pt", "GenMET_phi"]

import numpy as np
import operator

# make predictions to be used with keras.predict
def predict_inputs( event, sample, jet_lstm = False):
    flat_variables = np.array([[getattr( event, mva_variable) for mva_variable, _ in mva_variables]])
    if jet_lstm:
        lstm_jets_maxN = 10 #remove after retraining
        jet_vector_var = mva_vector_variables["mva_Jet"]
        jets = mva_vector_variables["mva_Jet"]["func"](event,sample=None)
        jets =  [ [ operator.itemgetter(varname)(jet) for varname in lstm_jetVarNames] for jet in jets[:lstm_jets_maxN] ]
        # zero padding
        jets += [ [0.]*len(lstm_jetVarNames)]*(max(0, lstm_jets_maxN-len(jets)))
        jets = np.array([jets])

        return [ flat_variables, jets ]
    else:
        return   flat_variables


# Dependence on TMB just for the sake of the example!

#define training samples for multiclassification
from TMB.Samples.nanoTuples_RunII_nanoAODv6_dilep_pp import * 
#use only Summer16
training_samples = [ Summer16.TTZ, Summer16.DY]#, Summer16.TTW]  

assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"

# training selection

from TMB.Tools.cutInterpreter import cutInterpreter
selectionString = cutInterpreter.cutString( 'dilepVL' )
