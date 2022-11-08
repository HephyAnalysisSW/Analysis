import ROOT
import os
from math import sqrt

from Analysis.Tools.helpers import getObjFromFile
from Analysis.Tools.u_float import u_float


maps_ele =           {"medium": "egammaEffi_EGM2D_Medium.root",
                      "tight":  "egammaEffi_EGM2D_Tight.root",
                      "loose":  "egammaEffi_EGM2D_Loose.root",
                      "Vloose": "egammaEffi_EGM2D_VLoose.root",
                     }

maps_mu =            {"tight":  "NUM_LeptonMvaTight_DEN_TrackerMuons_abseta_pt.root",
                      "medium": "NUM_LeptonMvaMedium_DEN_TrackerMuons_abseta_pt.root",
                      "loose":  "NUM_LeptonMvaLoose_DEN_TrackerMuons_abseta_pt.root",
                      "Vloose": "NUM_LeptonMvaVLoose_DEN_TrackerMuons_abseta_pt.root",
                     }

class leptonSF:
    def __init__(self, year, ID = None):

        if year not in ["2016_HIPM", "2016", "2017", "2018" ]:
            raise Exception("Lepton SF for year %i not known"%year)

        self.dataDir = "$CMSSW_BASE/src/Analysis/Tools/data/leptonSFData/LeptonMva_v1"
        self.year = year

        self.SFmaps = {
            "elec" : {
                "SF" :  getObjFromFile(self.dataDir+"/"+self.year+"/"+maps_el[ID],"EGamma_SF2D"),
                "syst": getObjFromFile(self.dataDir+"/"+self.year+"/"+maps_el[ID],"sys"),
                "stat": getObjFromFile(self.dataDir+"/"+self.year+"/"+maps_el[ID],"stat"),
            },
            "muon" : {
                "SF" :  getObjFromFile(self.dataDir+"/"+self.year+"/"+maps_mu[ID],"NUM_LeptonMvaTight_DEN_TrackerMuons"),
                "syst": getObjFromFile(self.dataDir+"/"+self.year+"/"+maps_mu[ID],"syst"),
                "stat": getObjFromFile(self.dataDir+"/"+self.year+"/"+maps_mu[ID],"stat"),
            }
        }


    def getSF(self, pdgId, pt, eta, unc='syst', sigma=0):
        uncert = "syst"
        if unc == "stat":
            uncert = "stat"
        lepton = None
        if abs(pdgId)==11:
            lepton = "elec"
            if eta > 2.5:
                eta = 2.49
            if eta < -2.5:
                eta = -2.49
            if pt > 200:
                pt = 199

        elif abs(pdgId)==13:
            lepton = "muon"
            eta = abs(eta)
            if eta > 2.4:
                eta = 2.39
            if pt > 120:
                pt = 119
        else:
          raise Exception("Lepton SF for PdgId %i not known"%pdgId)

        etabin = self.SFmaps[lepton]["SF"].GetXaxis().FindBin(eta)
        ptbin  = self.SFmaps[lepton]["SF"].GetYaxis().FindBin(pt)
        SF = self.SFmaps[lepton]["SF"].GetBinContent(etabin, ptbin)
        err = self.SFmaps[lepton][uncert].GetBinContent(etabin, ptbin)
        return SF+sigma*err
