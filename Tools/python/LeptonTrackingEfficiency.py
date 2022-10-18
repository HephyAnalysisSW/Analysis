import ROOT
import os, math
from Analysis.Tools.u_float import *
from Analysis.Tools.helpers import getObjFromFile

# Logging
import logging
logger = logging.getLogger(__name__)

class LeptonTrackingEfficiency:
    def __init__(self, year=2016):

        if year not in [ 2016, 2017, 2018 ]:
            raise Exception("Lepton Reconstruction Eff for year %i not known"%year)

        self.dataDir = "$CMSSW_BASE/src/Analysis/Tools/data/leptonSFData/"
        self.year = year

        if self.year == 2016:
            ## EGamma POG, https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#Cut%20based%20photon%20identification
            ## 2016 Legacy re-reco recommandation - Scale factors for reccomended ID for 2016 (renamed root file)
            e_file          = 'e2016_EGamma_Run2016BtoH_passingRECO_Legacy2016.root'
            e_file_lowEt    = 'e2016_EGamma_Run2016BtoH_passingRECO_lowEt_Legacy2016.root'
            e_key           = "EGamma_SF2D"

        elif self.year == 2017:
            ## EGamma POG, https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#Cut%20based%20photon%20identification
            e_file          = 'e2017_egammaEffi_EGM2D.root'
            e_file_lowEt    = 'e2017_egammaEffi_EGM2D_low.root'
            e_key           = "EGamma_SF2D"
            
        elif self.year == 2018:
            ## EGamma POG, https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#Cut%20based%20photon%20identification
            e_file          = 'e2018_egammaEffi_EGM2D.root'
            e_file_lowEt    = 'e2018_egammaEffi_EGM2D.root'
            e_key           = "EGamma_SF2D"
            
        self.e_sf       = getObjFromFile( os.path.expandvars( os.path.join( self.dataDir, e_file       ) ), e_key )
        self.e_sf_lowEt = getObjFromFile( os.path.expandvars( os.path.join( self.dataDir, e_file_lowEt ) ), e_key )

        assert self.e_sf,       "Could not load ele SF histo %s from file %s."%( e_key, e_file )
        assert self.e_sf_lowEt, "Could not load ele SF histo %s from file %s."%( e_key, e_file_lowEt )

        self.e_ptMax       = self.e_sf.GetYaxis().GetXmax()
        self.e_ptMin       = self.e_sf.GetYaxis().GetXmin()
        self.e_ptMin_lowEt = self.e_sf_lowEt.GetYaxis().GetXmin()

        self.e_etaMax      = self.e_sf.GetXaxis().GetXmax()
        self.e_etaMin      = self.e_sf.GetXaxis().GetXmin()

        ## Muons
        # SFs are 1. https://hypernews.cern.ch/HyperNews/CMS/get/muon/1425/1.html
        
    def getSF(self, pdgId, pt, eta, sigma=0):

        if abs(pdgId) == 11:
            if eta >= self.e_etaMax:
                logger.warning( "Supercluster eta out of bounds: %3.2f (need %3.2f <= eta <=% 3.2f)", eta, self.e_etaMin, self.e_etaMax )
                eta = self.e_etaMax - 0.01
            if eta <= self.e_etaMin:
                logger.warning( "Supercluster eta out of bounds: %3.2f (need %3.2f <= eta <=% 3.2f)", eta, self.e_etaMin, self.e_etaMax )
                eta = self.e_etaMin + 0.01

            # this is a bit awkward because of the seperate SFs for low pt electrons
            if   pt >= self.e_ptMax:       pt = self.e_ptMax - 1 
            elif pt <= self.e_ptMin_lowEt: pt = self.e_ptMin_lowEt + 1

            if pt <= self.e_ptMin: sf_hist = self.e_sf_lowEt
            else:                  sf_hist = self.e_sf

            val    = sf_hist.GetBinContent( sf_hist.FindBin(eta, pt) )
            valErr = sf_hist.GetBinError(   sf_hist.FindBin(eta, pt) )
            
            return val + sigma*valErr

        elif abs(pdgId) == 13:
            return 1

        else:
            raise ValueError( "Lepton pdgId %i neither electron or muon"%pdgId )



if __name__ == "__main__":

    sigma = 0
    print("2016")
    LSF = LeptonTrackingEfficiency(year=2016)
    print(LSF.getSF(11, 10, 1, sigma=sigma))
    print(LSF.getSF(11, 10, -1, sigma=sigma))
    print(LSF.getSF(13, 10, 1, sigma=sigma))
    print(LSF.getSF(13, 10, -1, sigma=sigma))

    print(LSF.getSF(11, 200, 1, sigma=sigma))
    print(LSF.getSF(11, 200, -1, sigma=sigma))
    print(LSF.getSF(13, 200, 1, sigma=sigma))
    print(LSF.getSF(13, 200, -1, sigma=sigma))

    print(LSF.getSF(11, 10, 2.5, sigma=sigma))
    print(LSF.getSF(11, 10, -2.5, sigma=sigma))
    print(LSF.getSF(13, 10, 2.5, sigma=sigma))
    print(LSF.getSF(13, 10, -2.5, sigma=sigma))

    print(LSF.getSF(11, 200, 2.5, sigma=sigma))
    print(LSF.getSF(11, 200, -2.5, sigma=sigma))
    print(LSF.getSF(13, 200, 2.5, sigma=sigma))
    print(LSF.getSF(13, 200, -2.5, sigma=sigma))

    print("2017")
    LSF = LeptonTrackingEfficiency(year=2017)
    print(LSF.getSF(11, 10, 1, sigma=sigma))
    print(LSF.getSF(11, 10, -1, sigma=sigma))
    print(LSF.getSF(13, 10, 1, sigma=sigma))
    print(LSF.getSF(13, 10, -1, sigma=sigma))

    print(LSF.getSF(11, 200, 1, sigma=sigma))
    print(LSF.getSF(11, 200, -1, sigma=sigma))
    print(LSF.getSF(13, 200, 1, sigma=sigma))
    print(LSF.getSF(13, 200, -1, sigma=sigma))

    print(LSF.getSF(11, 10, 2.5, sigma=sigma))
    print(LSF.getSF(11, 10, -2.5, sigma=sigma))
    print(LSF.getSF(13, 10, 2.5, sigma=sigma))
    print(LSF.getSF(13, 10, -2.5, sigma=sigma))

    print(LSF.getSF(11, 200, 2.5, sigma=sigma))
    print(LSF.getSF(11, 200, -2.5, sigma=sigma))
    print(LSF.getSF(13, 200, 2.5, sigma=sigma))
    print(LSF.getSF(13, 200, -2.5, sigma=sigma))

    print("2018")
    LSF = LeptonTrackingEfficiency(year=2018)
    print(LSF.getSF(11, 10, 1, sigma=sigma))
    print(LSF.getSF(11, 10, -1, sigma=sigma))
    print(LSF.getSF(13, 10, 1, sigma=sigma))
    print(LSF.getSF(13, 10, -1, sigma=sigma))

    print(LSF.getSF(11, 200, 1, sigma=sigma))
    print(LSF.getSF(11, 200, -1, sigma=sigma))
    print(LSF.getSF(13, 200, 1, sigma=sigma))
    print(LSF.getSF(13, 200, -1, sigma=sigma))

    print(LSF.getSF(11, 10, 2.5, sigma=sigma))
    print(LSF.getSF(11, 10, -2.5, sigma=sigma))
    print(LSF.getSF(13, 10, 2.5, sigma=sigma))
    print(LSF.getSF(13, 10, -2.5, sigma=sigma))

    print(LSF.getSF(11, 200, 2.5, sigma=sigma))
    print(LSF.getSF(11, 200, -2.5, sigma=sigma))
    print(LSF.getSF(13, 200, 2.5, sigma=sigma))
    print(LSF.getSF(13, 200, -2.5, sigma=sigma))

