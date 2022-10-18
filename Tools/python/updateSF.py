import os

commands = []

#Leptons
#leptonSF_datapath = os.path.expandvars('$CMSSW_BASE/src/Analysis/Tools/data/leptonSFData/')
#photonSF_datapath = os.path.expandvars('$CMSSW_BASE/src/Analysis/Tools/data/photonSFData/')
leptonSF_datapath = os.path.expandvars('$PWD/../data/leptonSFData/')
photonSF_datapath = os.path.expandvars('$PWD/../data/photonSFData/')

#Muons 2016
muonID_16_BCDEF_git  = 'https://gitlab.cern.ch/cms-muonPOG/MuonReferenceEfficiencies/raw/master/EfficienciesStudies/2016_legacy_rereco/systematic/RunBCDEF_SF_ID.root'
muonID_16_BCDEF_file = os.path.join( leptonSF_datapath, 'muon2016_RunBCDEF_SF_ID.root' )
commands.append( 'wget -O %s "%s"'%(muonID_16_BCDEF_file, muonID_16_BCDEF_git) )

muonISO_16_BCDEF_git  = 'https://gitlab.cern.ch/cms-muonPOG/MuonReferenceEfficiencies/raw/master/EfficienciesStudies/2016_legacy_rereco/systematic/RunBCDEF_SF_ISO.root'
muonISO_16_BCDEF_file = os.path.join( leptonSF_datapath, 'muon2016_RunBCDEF_SF_ISO.root' )
commands.append( 'wget -O %s "%s"'%(muonISO_16_BCDEF_file, muonISO_16_BCDEF_git) )

muonID_16_GH_git  = 'https://gitlab.cern.ch/cms-muonPOG/MuonReferenceEfficiencies/raw/master/EfficienciesStudies/2016_legacy_rereco/systematic/RunGH_SF_ID.root'
muonID_16_GH_file = os.path.join( leptonSF_datapath, 'muon2016_RunGH_SF_ID.root' )
commands.append( 'wget -O %s "%s"'%(muonID_16_GH_file, muonID_16_GH_git) )

muonISO_16_GH_git  = 'https://gitlab.cern.ch/cms-muonPOG/MuonReferenceEfficiencies/raw/master/EfficienciesStudies/2016_legacy_rereco/systematic/RunGH_SF_ISO.root'
muonISO_16_GH_file = os.path.join( leptonSF_datapath, 'muon2016_RunGH_SF_ISO.root' )
commands.append( 'wget -O %s "%s"'%(muonISO_16_GH_file, muonISO_16_GH_git) )

#Muons 2017
muonID_17_BCDEF_git  = 'https://twiki.cern.ch/twiki/pub/CMS/MuonReferenceEffs2017/RunBCDEF_SF_ID_syst.root'
muonID_17_BCDEF_file = os.path.join( leptonSF_datapath, 'muon2017_RunBCDEF_SF_ID.root' )
commands.append( 'wget -O %s "%s"'%(muonID_17_BCDEF_file, muonID_17_BCDEF_git) )

muonISO_17_BCDEF_git  = 'https://twiki.cern.ch/twiki/pub/CMS/MuonReferenceEffs2017/RunBCDEF_SF_ISO_syst.root'
muonISO_17_BCDEF_file = os.path.join( leptonSF_datapath, 'muon2017_RunBCDEF_SF_ISO.root' )
commands.append( 'wget -O %s "%s"'%(muonISO_17_BCDEF_file, muonISO_17_BCDEF_git) )

#Muons 2018
muonID_18_ABCD_git  = 'https://gitlab.cern.ch/cms-muonPOG/MuonReferenceEfficiencies/raw/master/EfficienciesStudies/2018/rootfiles/RunABCD_SF_ID.root'
muonID_18_ABCD_file = os.path.join( leptonSF_datapath, 'mu2018_RunABCD_SF_ID.root' )
commands.append( 'wget -O %s "%s"'%(muonID_18_ABCD_file, muonID_18_ABCD_git) )

muonISO_18_ABCD_git  = 'https://gitlab.cern.ch/cms-muonPOG/MuonReferenceEfficiencies/raw/master/EfficienciesStudies/2018/rootfiles/RunABCD_SF_ISO.root'
muonISO_18_ABCD_file = os.path.join( leptonSF_datapath, 'mu2018_RunABCD_SF_ISO.root' )
commands.append( 'wget -O %s "%s"'%(muonISO_18_ABCD_file, muonISO_18_ABCD_git) )


#Electrons 2016
eID_16_med_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/2016LegacyReReco_ElectronMedium.root'
eID_16_med_file = os.path.join( leptonSF_datapath, 'e2016_LegacyReReco_ElectronMedium.root' )
commands.append( 'wget -O %s "%s"'%(eID_16_med_file, eID_16_med_git) )

eID_16_tight_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/2016LegacyReReco_ElectronTight.root'
eID_16_tight_file = os.path.join( leptonSF_datapath, 'e2016_LegacyReReco_ElectronTight.root' )
commands.append( 'wget -O %s "%s"'%(eID_16_tight_file, eID_16_tight_git) )

eRec_16_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root'
eRec_16_file = os.path.join( leptonSF_datapath, 'e2016_EGamma_Run2016BtoH_passingRECO_Legacy2016.root' )
commands.append( 'wget -O %s "%s"'%(eRec_16_file, eRec_16_git) )

eRec_16_lowEt_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/EGM2D_BtoH_low_RecoSF_Legacy2016.root'
eRec_16_lowEt_file = os.path.join( leptonSF_datapath, 'e2016_EGamma_Run2016BtoH_passingRECO_lowEt_Legacy2016.root' )
commands.append( 'wget -O %s "%s"'%(eRec_16_lowEt_file, eRec_16_lowEt_git) )

#Electrons 2017
eID_17_med_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/2017_ElectronMedium.root'
eID_17_med_file = os.path.join( leptonSF_datapath, 'e2017_ElectronMediumCutBased.root' )
commands.append( 'wget -O %s "%s"'%(eID_17_med_file, eID_17_med_git) )

eID_17_tight_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/2017_ElectronTight.root'
eID_17_tight_file = os.path.join( leptonSF_datapath, 'e2017_ElectronTight.root' )
commands.append( 'wget -O %s "%s"'%(eID_17_tight_file, eID_17_tight_git) )

eRec_17_git  = 'https://twiki.cern.ch/twiki/pub/CMS/Egamma2017DataRecommendations/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root'
eRec_17_file = os.path.join( leptonSF_datapath, 'e2017_egammaEffi_EGM2D.root' )
commands.append( 'wget -O %s "%s"'%(eRec_17_file, eRec_17_git) )

eRec_17_lowEt_git  = 'https://twiki.cern.ch/twiki/pub/CMS/Egamma2017DataRecommendations/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root'
eRec_17_lowEt_file = os.path.join( leptonSF_datapath, 'e2017_egammaEffi_EGM2D_low.root' )
commands.append( 'wget -O %s "%s"'%(eRec_17_lowEt_file, eRec_17_lowEt_git) )

#Electrons 2018
eID_18_med_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/2018_ElectronMedium.root'
eID_18_med_file = os.path.join( leptonSF_datapath, 'e2018_ElectronMedium.root' )
commands.append( 'wget -O %s "%s"'%(eID_18_med_file, eID_18_med_git) )

eID_18_tight_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/2018_ElectronTight.root'
eID_18_tight_file = os.path.join( leptonSF_datapath, 'e2018_ElectronTight.root' )
commands.append( 'wget -O %s "%s"'%(eID_18_tight_file, eID_18_tight_git) )

# change to ONE file for pt>10 in leptonTracking!!!
eRec_18_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/egammaEffi.txt_EGM2D_updatedAll.root'
eRec_18_file = os.path.join( leptonSF_datapath, 'e2018_egammaEffi_EGM2D.root' )
commands.append( 'wget -O %s "%s"'%(eRec_18_file, eRec_18_git) )

eRec_18_lowEt_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/egammaEffi.txt_EGM2D_updatedAll.root'
eRec_18_lowEt_file = os.path.join( leptonSF_datapath, 'e2018_egammaEffi_EGM2D_low.root' )
commands.append( 'wget -O %s "%s"'%(eRec_18_lowEt_file, eRec_18_lowEt_git) )






#Photons 2016
gID_16_med_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/80X_2016_Medium_photons.root'
gID_16_med_file = os.path.join( photonSF_datapath, 'g2016_LegacyReReco_PhotonCutBasedMedium.root' )
commands.append( 'wget -O %s "%s"'%(gID_16_med_file, gID_16_med_git) )

#gRec_16_git  = ''
#gRec_16_file = os.path.join( photonSF_datapath, 'g2016_EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root' )
#commands.append( 'wget -O %s "%s"'%(gRec_16_file, gRec_16_git) )

gVeto_16_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/ScalingFactors_80X_Summer16.root'
gVeto_16_file = os.path.join( photonSF_datapath, 'g2016_ScalingFactors_80X_Summer16.root' )
commands.append( 'wget -O %s "%s"'%(gVeto_16_file, gVeto_16_git) )

#Photons 2017
gID_17_med_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/2017_PhotonsMedium.root'
gID_17_med_file = os.path.join( photonSF_datapath, 'g2017_PhotonsMedium.root' )
commands.append( 'wget -O %s "%s"'%(gID_17_med_file, gID_17_med_git) )

#gRec_17_git  = ''
#gRec_17_file = os.path.join( photonSF_datapath, '' )
#commands.append( 'wget -O %s "%s"'%(gRec_17_file, gRec_17_git) )

#gVeto_17_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/PixelSeed_ScaleFactors_2017.root'
gVeto_17_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/CSEV_ScaleFactors_2017.root'
gVeto_17_file = os.path.join( photonSF_datapath, 'g2017_PixelSeedVeto_ScaleFactors.root' )
commands.append( 'wget -O %s "%s"'%(gVeto_17_file, gVeto_17_git) )

#Photons 2018
gID_18_med_git  = 'https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/2018_PhotonsMedium.root'
gID_18_med_file = os.path.join( photonSF_datapath, 'g2018_PhotonsMedium.root' )
commands.append( 'wget -O %s "%s"'%(gID_18_med_file, gID_18_med_git) )

#gRec_18_git  = ''
#gRec_18_file = os.path.join( photonSF_datapath, '' )
#commands.append( 'wget -O %s "%s"'%(gRec_18_file, gRec_18_git) )

#gVeto_18_git  = ''
#gVeto_18_file = os.path.join( photonSF_datapath, 'g2018_PixelSeedVeto_ScaleFactors.root' )
#commands.append( 'wget -O %s "%s"'%(gVeto_18_file, gVeto_18_git) )


for command in commands:
#    os.system(command)
    print(command)
