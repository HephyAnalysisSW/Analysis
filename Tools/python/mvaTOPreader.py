import os
import xgboost as xgb
import numpy as np
from math import log, log10
import logging
logger = logging.getLogger(__name__)

class mvaTOPreader:
    def __init__(self, year):
        yearstring = 'UL18'
        if year == "UL2016":
            yearstring = 'UL16'
        elif year == "UL2016_preVFP": 
            yearstring = 'UL16APV'
        elif year == "UL2017":
            yearstring = 'UL17'
        elif year == "UL2018":
            yearstring = 'UL18'
            
        cmsswbase=os.environ["CMSSW_BASE"]
        directory = cmsswbase+"/src/tWZ/Tools/data/mvaWeights/"
        
        # Load Electron weights
        # v1
        self.bst_el = xgb.Booster() 
        self.bst_el.load_model(directory+'el_TOP'+yearstring+'_XGB.weights.bin') 
        self.bst_elv2 = xgb.Booster() 
        #v2
        self.bst_elv2.load_model(directory+'el_TOPv2'+yearstring+'_XGB.weights.bin') 
        
        # Load Muon weights
        # v1
        self.bst_mu = xgb.Booster()    
        self.bst_mu.load_model(directory+'mu_TOP'+yearstring+'_XGB.weights.bin') 
        # v2
        self.bst_muv2 = xgb.Booster()    
        self.bst_muv2.load_model(directory+'mu_TOPv2'+yearstring+'_XGB.weights.bin')         
        
        # Working points
        self.WPs    = [0.20, 0.41, 0.64, 0.81]
        self.WPs_v2 = [0.59, 0.81, 0.90, 0.94]

        

    def getmvaTOPScore(self, lep):
        mvaScore = -1
        mvaScorev2 = -1
        if abs(lep['pdgId']) == 11:
            features = np.array([[
                lep['pt'], 
                lep['eta'], 
                ord(lep['jetNDauCharged']), # jetNDauChargedMVASel
                lep['miniPFRelIso_chg'], # miniRelIsoCharged
                lep['miniPFRelIso_all']-lep['miniPFRelIso_chg'], # miniRelIsoNeutralVanilla
                lep['jetPtRelv2'],
                lep['jetPtRatio'], # jetPtRatioVanilla 
                lep['pfRelIso03_all'], # relIso0p3Vanilla
                lep['jetBTag'],
            	lep['sip3d'],
            	log(abs(lep['dxy'])),
            	log(abs(lep['dz'])),
                lep['mvaFall17V2noIso'], # eleMvaFall17v2
                ord(lep['lostHits']), # eleMissingHits
            ]])
            dtest = xgb.DMatrix(features)
            mvaScore = self.bst_el.predict(dtest)
            mvaScorev2 = self.bst_elv2.predict(dtest)
        elif abs(lep['pdgId']) == 13:
            features = np.array([[
                lep['pt'], 
                lep['eta'], 
                ord(lep['jetNDauCharged']), # jetNDauChargedMVASel
                lep['miniPFRelIso_chg'], # miniRelIsoCharged
                lep['miniPFRelIso_all']-lep['miniPFRelIso_chg'], # miniRelIsoNeutralVanilla
                lep['jetPtRelv2'],
                lep['jetPtRatio'], # jetPtRatioVanilla
                lep['pfRelIso03_all'], # relIso0p3DBVanilla
                lep['jetBTag'],
            	lep['sip3d'],
            	log(abs(lep['dxy'])), 
            	log(abs(lep['dz'])), 
                lep['segmentComp'], # segComp
            ]])
            dtest = xgb.DMatrix(features)
            mvaScore = self.bst_mu.predict(dtest)
            mvaScorev2 = self.bst_muv2.predict(dtest)
        
        # Now get Working point for v1 and v2
        WPv1 = 0
        for wp in self.WPs:
            if mvaScore > wp:
                WPv1 += 1
        WPv2 = 0
        for wp in self.WPs_v2:
            if mvaScorev2 > wp:
                WPv2 += 1
        ####
        
        return mvaScore, WPv1, mvaScorev2, WPv2
