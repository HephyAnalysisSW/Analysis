import ROOT
import os
from math import sqrt

from Analysis.Tools.helpers import getObjFromFile
from Analysis.Tools.u_float import u_float

# 2016 Lumi Ratios
lumiRatio2016_BCDEF = 19.695422959 / 35.921875595
lumiRatio2016_GH    = 16.226452636 / 35.921875595

keys_mu2016_BCDEF = { 
                      "medium":[( "mu2016_RunBCDEF_SF_ID.root",  "NUM_MediumID_DEN_genTracks_eta_pt"   ),
                                ( "mu2016_RunBCDEF_SF_ISO.root", "NUM_TightRelIso_DEN_MediumID_eta_pt" )],
                      "medium_stat":[( "mu2016_RunBCDEF_SF_ID.root",  "NUM_MediumID_DEN_genTracks_eta_pt_stat"   ),
                                ( "mu2016_RunBCDEF_SF_ISO.root", "NUM_TightRelIso_DEN_MediumID_eta_pt_stat" )],
                      "medium_syst":[( "mu2016_RunBCDEF_SF_ID.root",  "NUM_MediumID_DEN_genTracks_eta_pt_syst"   ),
                                ( "mu2016_RunBCDEF_SF_ISO.root", "NUM_TightRelIso_DEN_MediumID_eta_pt_syst" )],
                      "tight": [( "mu2016_RunBCDEF_SF_ID.root",  "NUM_TightID_DEN_genTracks_eta_pt"   ),
                                ( "mu2016_RunBCDEF_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt" )],
                      "tight_stat": [( "mu2016_RunBCDEF_SF_ID.root",  "NUM_TightID_DEN_genTracks_eta_pt_stat"   ),
                                     ( "mu2016_RunBCDEF_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt_stat" )],
                      "tight_syst": [( "mu2016_RunBCDEF_SF_ID.root",  "NUM_TightID_DEN_genTracks_eta_pt_syst"   ),
                                     ( "mu2016_RunBCDEF_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt_syst" )],
                    }

keys_mu2016_GH    = { 
                      "medium":[( "mu2016_RunGH_SF_ID.root",     "NUM_MediumID_DEN_genTracks_eta_pt"   ),
                                ( "mu2016_RunGH_SF_ISO.root",    "NUM_TightRelIso_DEN_MediumID_eta_pt" )],
                      "medium_stat":[( "mu2016_RunGH_SF_ID.root",     "NUM_MediumID_DEN_genTracks_eta_pt_stat"   ),
                                ( "mu2016_RunGH_SF_ISO.root",    "NUM_TightRelIso_DEN_MediumID_eta_pt_stat" )],
                      "medium_syst":[( "mu2016_RunGH_SF_ID.root",     "NUM_MediumID_DEN_genTracks_eta_pt_syst"   ),
                                ( "mu2016_RunGH_SF_ISO.root",    "NUM_TightRelIso_DEN_MediumID_eta_pt_syst" )],
                      "tight": [( "mu2016_RunGH_SF_ID.root",     "NUM_TightID_DEN_genTracks_eta_pt"   ),
                                ( "mu2016_RunGH_SF_ISO.root",    "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt" )],
                      "tight_stat": [( "mu2016_RunGH_SF_ID.root",     "NUM_TightID_DEN_genTracks_eta_pt_stat"   ),
                                     ( "mu2016_RunGH_SF_ISO.root",    "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt_stat" )],
                      "tight_syst": [( "mu2016_RunGH_SF_ID.root",     "NUM_TightID_DEN_genTracks_eta_pt_syst"   ),
                                     ( "mu2016_RunGH_SF_ISO.root",    "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt_syst" )],
                    }

keys_ele2016      = { "medium":[( "e2016_LegacyReReco_ElectronMedium_Fall17V2.root", "EGamma_SF2D" )],
                      "tight": [( "e2016_LegacyReReco_ElectronTight_Fall17V2.root", "EGamma_SF2D" )],
                    }

keys_mu2017       = { 
                      "medium":[( "mu2017_RunBCDEF_SF_ID_syst.root",  "NUM_MediumID_DEN_genTracks_pt_abseta"   ),
                                ( "mu2017_RunBCDEF_SF_ISO_syst.root", "NUM_TightRelIso_DEN_MediumID_pt_abseta" )],
                      "medium_stat":[( "mu2017_RunBCDEF_SF_ID_syst.root",  "NUM_MediumID_DEN_genTracks_pt_abseta_stat"   ),
                                ( "mu2017_RunBCDEF_SF_ISO_syst.root", "NUM_TightRelIso_DEN_MediumID_pt_abseta_stat" )],
                      "medium_syst":[( "mu2017_RunBCDEF_SF_ID_syst.root",  "NUM_MediumID_DEN_genTracks_pt_abseta_syst"   ),
                                ( "mu2017_RunBCDEF_SF_ISO_syst.root", "NUM_TightRelIso_DEN_MediumID_pt_abseta_syst" )],
                      "tight": [( "mu2017_RunBCDEF_SF_ID_syst.root",  "NUM_TightID_DEN_genTracks_pt_abseta"   ),
                                ( "mu2017_RunBCDEF_SF_ISO_syst.root", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta" )],
                      "tight_stat": [( "mu2017_RunBCDEF_SF_ID_syst.root",  "NUM_TightID_DEN_genTracks_pt_abseta_stat"   ),
                                     ( "mu2017_RunBCDEF_SF_ISO_syst.root", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_stat" )],
                      "tight_syst": [( "mu2017_RunBCDEF_SF_ID_syst.root",  "NUM_TightID_DEN_genTracks_pt_abseta_syst"   ),
                                     ( "mu2017_RunBCDEF_SF_ISO_syst.root", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst" )],
                    }

keys_ele2017      = { "medium":[( "e2017_ElectronMediumCutBased.root", "EGamma_SF2D" )],
                      "tight": [( "e2017_ElectronTight.root", "EGamma_SF2D" )],
                    }

keys_mu2018       = { 
                      "medium":[( "mu2018_RunABCD_SF_ID.root",  "NUM_MediumID_DEN_TrackerMuons_pt_abseta"   ),
                                ( "mu2018_RunABCD_SF_ISO.root", "NUM_TightRelIso_DEN_MediumID_pt_abseta" )],
                      "medium_stat":[( "mu2018_RunABCD_SF_ID.root",  "NUM_MediumID_DEN_TrackerMuons_pt_abseta_stat"   ),
                                ( "mu2018_RunABCD_SF_ISO.root", "NUM_TightRelIso_DEN_MediumID_pt_abseta_stat" )],
                      "medium_syst":[( "mu2018_RunABCD_SF_ID.root",  "NUM_MediumID_DEN_TrackerMuons_pt_abseta_syst"   ),
                                ( "mu2018_RunABCD_SF_ISO.root", "NUM_TightRelIso_DEN_MediumID_pt_abseta_syst" )],
                      "tight": [( "mu2018_RunABCD_SF_ID.root",  "NUM_TightID_DEN_TrackerMuons_pt_abseta"   ),
                                ( "mu2018_RunABCD_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta" )],
                      "tight_stat": [( "mu2018_RunABCD_SF_ID.root",  "NUM_TightID_DEN_TrackerMuons_pt_abseta_stat"   ),
                                     ( "mu2018_RunABCD_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_stat" )],
                      "tight_syst": [( "mu2018_RunABCD_SF_ID.root",  "NUM_TightID_DEN_TrackerMuons_pt_abseta_syst"   ),
                                     ( "mu2018_RunABCD_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst" )],
                    }


keys_ele2018      = { "medium":[( "e2018_ElectronMedium.root", "EGamma_SF2D" )],
                      "tight": [( "e2018_ElectronTight.root", "EGamma_SF2D" )],
                    }


class LeptonSF:

    def __init__(self, year=2016, ID=None):

        if year not in [ 2016, 2017, 2018 ]:
            raise Exception("Lepton SF for year %i not known"%year)

        self.dataDir = "$CMSSW_BASE/src/Analysis/Tools/data/leptonSFData"
        self.year    = year

        if year == 2016:

            if not ID in keys_mu2016_BCDEF.keys():
                raise Exception("Don't know ID %s"%ID)

            if not ID in keys_mu2016_GH.keys():
                raise Exception("Don't know ID %s"%ID)

            if not ID in keys_ele2016.keys():
                raise Exception("Don't know ID %s"%ID)

            self.mu_BCDEF = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2016_BCDEF[ID]              ]
            self.mu_GH    = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2016_GH[ID]                 ]
            self.mu_BCDEF_stat = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2016_BCDEF[ID+"_stat"] ]
            self.mu_GH_stat    = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2016_GH[ID+"_stat"]    ]
            self.mu_BCDEF_syst = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2016_BCDEF[ID+"_syst"] ]
            self.mu_GH_syst    = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2016_GH[ID+"_syst"]    ]
            self.ele      = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_ele2016[ID]                   ]

            for effMap in self.mu_BCDEF + self.mu_GH + self.ele + self.mu_BCDEF_stat + self.mu_BCDEF_syst + self.mu_GH_stat + self.mu_GH_syst: assert effMap

        elif year == 2017:

            if not ID in keys_mu2017.keys():
                raise Exception("Don't know ID %s"%ID)

            if not ID in keys_ele2017.keys():
                raise Exception("Don't know ID %s"%ID)

            self.mu       = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2017[ID]         ]
            self.mu_stat  = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2017[ID+"_stat"] ]
            self.mu_syst  = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2017[ID+"_syst"] ]
            self.ele      = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_ele2017[ID]        ]

            for effMap in self.mu + self.ele + self.mu_stat + self.mu_syst: assert effMap

        elif year == 2018:

            if not ID in keys_mu2018.keys():
                raise Exception("Don't know ID %s"%ID)

            if not ID in keys_ele2018.keys():
                raise Exception("Don't know ID %s"%ID)

            self.mu       = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2018[ID]         ]
            self.mu_stat  = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2018[ID+"_stat"] ]
            self.mu_syst  = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_mu2018[ID+"_syst"] ]
            self.ele      = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in keys_ele2018[ID]        ]

            for effMap in self.ele + self.mu + self.mu_stat + self.mu_syst: assert effMap

    def getPartialSF( self, effMap, pt, eta, reversed=False ):
        x = eta if not reversed else pt
        y = pt  if not reversed else eta
        sf  = effMap.GetBinContent( effMap.FindBin(x, y) )
        err = effMap.GetBinError(   effMap.FindBin(x, y) )
        return u_float(sf, err)

    def mult( self, list ):
        res = list[0]
        for i in list[1:]: res = res*i
        return res

    def getSF(self, pdgId, pt, eta, sigma=0, unc="nominal"):

        if abs(pdgId) not in [11,13]:
            raise Exception("Lepton SF for PdgId %i not known"%pdgId)

        if not unc in ["nominal", "stat", "syst"]:
            raise Exception("Don't know uncertainty %s"%unc)

        if abs(pdgId) == 11 and unc != "nominal":
            raise Exception("Stat and syst uncertainty only implemented for muons")

        if self.year == 2016 and abs(pdgId) == 13:
            if   pt  >=  120: pt  =   119
            if   pt  <=   20: pt  =    21
            if   eta >=  2.4: eta =  2.39 
            elif eta <= -2.4: eta = -2.39 

            if unc == "nominal":
                mu_BCDEF = self.mu_BCDEF
                mu_GH    = self.mu_GH
            elif unc == "stat":
                mu_BCDEF = self.mu_BCDEF_stat
                mu_GH    = self.mu_GH_stat
            elif unc == "syst":
                mu_BCDEF = self.mu_BCDEF_syst
                mu_GH    = self.mu_GH_syst

            sf_BCDEF = self.mult( [self.getPartialSF(effMap, pt, eta) for effMap in mu_BCDEF] )
            sf_GH    = self.mult( [self.getPartialSF(effMap, pt, eta) for effMap in mu_GH] )
            sf       = sf_BCDEF*lumiRatio2016_BCDEF + sf_GH*lumiRatio2016_GH
#            print "B", sf_BCDEF, "G", sf_GH, "SF", sf

#            print "BID", [self.getPartialSF(effMap, pt, eta) for effMap in self.mu_BCDEF][0]
#            print "BIso", [self.getPartialSF(effMap, pt, eta) for effMap in self.mu_BCDEF][1]

#            print "GID", [self.getPartialSF(effMap, pt, eta) for effMap in self.mu_GH][0]
#            print "GIso", [self.getPartialSF(effMap, pt, eta) for effMap in self.mu_GH][1]

        elif self.year == 2018 and abs(pdgId) == 13:
            absEta = abs(eta)
            if pt     >= 120: pt     = 119
            if pt     <= 20:  pt     = 21
            if absEta >= 2.4: absEta = 2.39 

            if unc == "nominal":
                mu = self.mu
            elif unc == "stat":
                mu = self.mu_stat
            elif unc == "syst":
                mu = self.mu_syst

            sf = self.mult( [ self.getPartialSF( effMap, pt, absEta, reversed=True ) for effMap in mu ] )

        else:
            if abs(pdgId) == 13:
                absEta = abs(eta)
                if pt     >= 120: pt     =  119
                if pt     <=  20: pt     =   21
                if absEta >= 2.4: absEta = 2.39 

                if unc == "nominal":
                    mu = self.mu
                elif unc == "stat":
                    mu = self.mu_stat
                elif unc == "syst":
                    mu = self.mu_syst

                sf = self.mult( [ self.getPartialSF( effMap, pt, absEta, reversed=True ) for effMap in mu ] )

            elif abs(pdgId) == 11:
                if   pt  >=  500: pt  =   499
                if   pt  <=   10: pt  =    11
                if   eta >=  2.5: eta =  2.49 
                elif eta <= -2.5: eta = -2.49 

                sf = self.mult( [ self.getPartialSF( effMap, pt, eta ) for effMap in self.ele ] )

        return (1+sf.sigma*sigma)*sf.val


if __name__ == "__main__":

    sigma = 1
    print "2016, medium"
    LSF = LeptonSF(year=2016, ID="medium")
    print LSF.getSF(11, 10, 1, sigma=sigma)
    print LSF.getSF(11, 10, -1, sigma=sigma)
    print LSF.getSF(13, 10, 1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, -1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, 1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, -1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, 1, sigma=sigma, unc="stat")
    print LSF.getSF(13, 10, -1, sigma=sigma, unc="stat")

    print LSF.getSF(11, 200, 1, sigma=sigma)
    print LSF.getSF(11, 200, -1, sigma=sigma)
    print LSF.getSF(13, 200, 1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, -1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, 1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, -1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, 1, sigma=sigma, unc="stat")
    print LSF.getSF(13, 200, -1, sigma=sigma, unc="stat")

    print LSF.getSF(11, 10, 2.5, sigma=sigma)
    print LSF.getSF(11, 10, -2.5, sigma=sigma)
    print LSF.getSF(13, 10, 2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, -2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, 2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, -2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, 2.5, sigma=sigma, unc="stat")
    print LSF.getSF(13, 10, -2.5, sigma=sigma, unc="stat")

    print LSF.getSF(11, 200, 2.5, sigma=sigma)
    print LSF.getSF(11, 200, -2.5, sigma=sigma)
    print LSF.getSF(13, 200, 2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, -2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, 2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, -2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, 2.5, sigma=sigma, unc="stat")
    print LSF.getSF(13, 200, -2.5, sigma=sigma, unc="stat")

    print "2017, medium"
    LSF = LeptonSF(year=2017, ID="medium")
    print LSF.getSF(11, 10, 1, sigma=sigma)
    print LSF.getSF(11, 10, -1, sigma=sigma)
    print LSF.getSF(13, 10, 1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, -1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, 1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, -1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, 1, sigma=sigma, unc="stat")
    print LSF.getSF(13, 10, -1, sigma=sigma, unc="stat")

    print LSF.getSF(11, 200, 1, sigma=sigma)
    print LSF.getSF(11, 200, -1, sigma=sigma)
    print LSF.getSF(13, 200, 1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, -1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, 1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, -1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, 1, sigma=sigma, unc="stat")
    print LSF.getSF(13, 200, -1, sigma=sigma, unc="stat")

    print LSF.getSF(11, 10, 2.5, sigma=sigma)
    print LSF.getSF(11, 10, -2.5, sigma=sigma)
    print LSF.getSF(13, 10, 2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, -2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, 2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, -2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, 2.5, sigma=sigma, unc="stat")
    print LSF.getSF(13, 10, -2.5, sigma=sigma, unc="stat")

    print LSF.getSF(11, 200, 2.5, sigma=sigma)
    print LSF.getSF(11, 200, -2.5, sigma=sigma)
    print LSF.getSF(13, 200, 2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, -2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, 2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, -2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, 2.5, sigma=sigma, unc="stat")
    print LSF.getSF(13, 200, -2.5, sigma=sigma, unc="stat")

    print "2018, medium"
    LSF = LeptonSF(year=2018, ID="medium")
    print LSF.getSF(11, 10, 1, sigma=sigma)
    print LSF.getSF(11, 10, -1, sigma=sigma)
    print LSF.getSF(13, 10, 1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, -1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, 1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, -1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, 1, sigma=sigma, unc="stat")
    print LSF.getSF(13, 10, -1, sigma=sigma, unc="stat")

    print LSF.getSF(11, 200, 1, sigma=sigma)
    print LSF.getSF(11, 200, -1, sigma=sigma)
    print LSF.getSF(13, 200, 1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, -1, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, 1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, -1, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, 1, sigma=sigma, unc="stat")
    print LSF.getSF(13, 200, -1, sigma=sigma, unc="stat")

    print LSF.getSF(11, 10, 2.5, sigma=sigma)
    print LSF.getSF(11, 10, -2.5, sigma=sigma)
    print LSF.getSF(13, 10, 2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, -2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 10, 2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, -2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 10, 2.5, sigma=sigma, unc="stat")
    print LSF.getSF(13, 10, -2.5, sigma=sigma, unc="stat")

    print LSF.getSF(11, 200, 2.5, sigma=sigma)
    print LSF.getSF(11, 200, -2.5, sigma=sigma)
    print LSF.getSF(13, 200, 2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, -2.5, sigma=sigma, unc="nominal")
    print LSF.getSF(13, 200, 2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, -2.5, sigma=sigma, unc="syst")
    print LSF.getSF(13, 200, 2.5, sigma=sigma, unc="stat")
    print LSF.getSF(13, 200, -2.5, sigma=sigma, unc="stat")
