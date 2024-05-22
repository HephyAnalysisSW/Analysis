import ROOT
import os


from Analysis.Tools.helpers import getObjFromFile

#https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
#eta has to be supercluster eta

map_el = {"UL2016_preVFP":{"below20pt":"egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root", "above20pt":"egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root"},
          "UL2016":{"below20pt":"egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root", "above20pt":"egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root"},
          "UL2017":{"below20pt":"egammaEffi_ptBelow20.txt_EGM2D_UL2017.root", "above20pt":"egammaEffi_ptAbove20.txt_EGM2D_UL2017.root"},
          "UL2018":{"below20pt":"egammaEffi_ptBelow20.txt_EGM2D_UL2018.root", "above20pt":"egammaEffi_ptAbove20.txt_EGM2D_UL2018.root"},
        }

class ElectronRecoSF:
    def __init__(self, era):

        if era not in ["UL2016_preVFP", "UL2016", "UL2017", "UL2018" ]:
            raise Exception("Electron reco efficiency for era %i not known"%era)

        self.dataDir = "$CMSSW_BASE/src/Analysis/Tools/data/electronRecoData"
        self.era = era


    def getSF(self, pt, eta, sigma=0):
        
        if pt>500 or pt<10:
            raise Exception("Electron pt out of range %s"%pt)
        if pt >=20 : ptKey = "above20pt"
        else: ptKey = "below20pt"

        effHist = getObjFromFile(self.dataDir+"/"+map_el[self.era][ptKey],"EGamma_SF2D")
         
        if eta > 2.5:
            eta = 2.49
        if eta < -2.5:
            eta = -2.49

        etabin = effHist.GetXaxis().FindBin(eta)
        ptbin  = effHist.GetYaxis().FindBin(pt)
        SF = effHist.GetBinContent(etabin, ptbin)
        err = effHist.GetBinError(etabin, ptbin)


        return SF+sigma*err

if __name__ == '__main__':

    sf = ElectronRecoSF("UL2016_preVFP")
    print sf.getSF(pt=20, eta= 0.5)
    print sf.getSF(pt=50, eta= 0.5, sigma=1)
    print sf.getSF(pt=15, eta= 0.5, sigma=-1)
