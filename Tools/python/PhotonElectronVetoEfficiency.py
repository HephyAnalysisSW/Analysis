import ROOT
import os

from Analysis.Tools.helpers import getObjFromFile
from Analysis.Tools.u_float import u_float

# Logging
import logging
logger = logging.getLogger(__name__)

g16_keys = [( "g2016_ScalingFactors_80X_Summer16.root", "Scaling_Factors_CSEV_R9 Inclusive"   ),
            ( "g2016_ScalingFactors_80X_Summer16.root", "Scaling_Factors_HasPix_R9 Inclusive" )]

g17_keys = [( "g2017_CSEV_ScaleFactors_2017.root",      "Medium_ID"   ),
            ( "g2017_PixelSeed_ScaleFactors_2017.root", "Medium_ID" )]

g18_keys = [( "g2018_HasPix_2018.root", "eleVeto_SF"   ),
            ( "g2018_CSEV_2018.root",   "eleVeto_SF" )]

g18_keys_unc = [( "g2018_HasPix_2018.root", "eleVeto_Unc"   ),
                ( "g2018_CSEV_2018.root",   "eleVeto_Unc" )]

class PhotonElectronVetoEfficiency:
    def __init__( self, year ):

        if year not in [ 2016, 2017, 2018 ]:
            raise Exception("Photon Veto SF for year %i not known"%year)

        self.year = year

        if self.year == 2016:
            g_keys = g16_keys
        elif self.year == 2017:
            g_keys = g17_keys
        elif self.year == 2018:
            g_keys     = g18_keys
            g_keys_unc = g18_keys_unc

        self.dataDir = "$CMSSW_BASE/src/Analysis/Tools/data/photonSFData"
        self.g_sf = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in g_keys ]
        if self.year == 2018:
            self.g_sf_unc = [ getObjFromFile(os.path.expandvars(os.path.join(self.dataDir, file)), key) for (file, key) in g_keys_unc ]
            for effMap in self.g_sf_unc: assert effMap

        for effMap in self.g_sf: assert effMap

    def getPartialSF( self, effMap, pt, eta, effMap_unc=None ):
        if self.year == 2016:
            bin = effMap.FindBin(eta,pt)
        elif self.year == 2017:
            bin = 1 if abs(eta) < 1.479 else 4
        elif self.year == 2018:
            bin = effMap.FindBin(pt,eta)

        sf  = effMap.GetBinContent( bin )

        if self.year == 2018:
            err = effMap_unc.GetBinContent( bin )
        else:
            err = effMap.GetBinError( bin )

        return u_float(sf, err)

    def mult( self, list ):
        res = list[0]
        for i in list[1:]: res = res*i
        return res

    def getSF( self, pt, eta, sigma=0 ):
        if pt       >= 200: pt  = 199
        if abs(eta) >= 2.5: eta = 2.49 
        if self.year == 2018:
            sf = self.mult( [ self.getPartialSF( effMap, pt, abs(eta), self.g_sf_unc[i] ) for i, effMap in enumerate(self.g_sf) ] )
        else:
            sf = self.mult( [ self.getPartialSF( effMap, pt, abs(eta) ) for effMap in self.g_sf ] )

        return (1+sf.sigma*sigma)*sf.val


if __name__ == "__main__":

    sigma = 0
    print "2016"
    LSF = PhotonElectronVetoEfficiency(year=2016)
    print LSF.getSF(20, 1, sigma=sigma)
    print LSF.getSF(20, -1, sigma=sigma)
    print LSF.getSF(20, 1, sigma=sigma)
    print LSF.getSF(20, -1, sigma=sigma)

    print LSF.getSF(200, 1, sigma=sigma)
    print LSF.getSF(200, -1, sigma=sigma)
    print LSF.getSF(200, 1, sigma=sigma)
    print LSF.getSF(200, -1, sigma=sigma)

    print LSF.getSF(20, 1.4, sigma=sigma)
    print LSF.getSF(20, -1.4, sigma=sigma)
    print LSF.getSF(20, 1.4, sigma=sigma)
    print LSF.getSF(20, -1.4, sigma=sigma)

    print LSF.getSF(200, 1.4, sigma=sigma)
    print LSF.getSF(200, -1.4, sigma=sigma)
    print LSF.getSF(200, 1.4, sigma=sigma)
    print LSF.getSF(200, -1.4, sigma=sigma)

    print "2017"
    LSF = PhotonElectronVetoEfficiency(year=2017)
    print LSF.getSF(20, 1, sigma=sigma)
    print LSF.getSF(20, -1, sigma=sigma)
    print LSF.getSF(20, 1, sigma=sigma)
    print LSF.getSF(20, -1, sigma=sigma)

    print LSF.getSF(200, 1, sigma=sigma)
    print LSF.getSF(200, -1, sigma=sigma)
    print LSF.getSF(200, 1, sigma=sigma)
    print LSF.getSF(200, -1, sigma=sigma)

    print LSF.getSF(20, 1.4, sigma=sigma)
    print LSF.getSF(20, -1.4, sigma=sigma)
    print LSF.getSF(20, 1.4, sigma=sigma)
    print LSF.getSF(20, -1.4, sigma=sigma)

    print LSF.getSF(200, 1.4, sigma=sigma)
    print LSF.getSF(200, -1.4, sigma=sigma)
    print LSF.getSF(200, 1.4, sigma=sigma)
    print LSF.getSF(200, -1.4, sigma=sigma)

    print "2018"
    LSF = PhotonElectronVetoEfficiency(year=2018)
    print LSF.getSF(20, 1, sigma=sigma)
    print LSF.getSF(20, -1, sigma=sigma)
    print LSF.getSF(20, 1, sigma=sigma)
    print LSF.getSF(20, -1, sigma=sigma)

    print LSF.getSF(200, 1, sigma=sigma)
    print LSF.getSF(200, -1, sigma=sigma)
    print LSF.getSF(200, 1, sigma=sigma)
    print LSF.getSF(200, -1, sigma=sigma)

    print LSF.getSF(20, 1.4, sigma=sigma)
    print LSF.getSF(20, -1.4, sigma=sigma)
    print LSF.getSF(20, 1.4, sigma=sigma)
    print LSF.getSF(20, -1.4, sigma=sigma)

    print LSF.getSF(200, 1.4, sigma=sigma)
    print LSF.getSF(200, -1.4, sigma=sigma)
    print LSF.getSF(200, 1.4, sigma=sigma)
    print LSF.getSF(200, -1.4, sigma=sigma)

