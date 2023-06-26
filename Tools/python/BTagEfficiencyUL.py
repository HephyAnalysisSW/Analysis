''' Implementation of b-tagging reweighting
'''

# Standard imports
import ROOT, pickle, itertools, os
from operator import mul
from correctionlib import _core

# Logging
import logging
logger = logging.getLogger(__name__)

#binning in pt and eta
ptBorders = [20, 30, 50, 70, 100, 140, 200, 300, 600, 1000]
ptBins    = [ [ptBorders[i], ptBorders[i+1]] for i in range(len(ptBorders)-1) ]
ptBins   += [ [ptBorders[-1], -1] ]

etaBins2016 = [[0,2.4]]
etaBins2017 = [[0,2.4]]
etaBins2018 = [[0,2.5]]

def toFlavourKey(pdgId):
    if abs(pdgId)==5: return 0
    if abs(pdgId)==4: return 1
    return 2

#Method 1ab
#UL Files

##All efficiencies are now with medium WP
effFile2016postVFPULDeepCSV = 'TTGJets_2016_2j_1l_DeepB_eta_v3.pkl'
effFile2016preVFPULDeepCSV = 'TTGJets_2016APV_2j_1l_DeepB_eta_v2.pkl'
effFile2017ULDeepCSV = 'TTGJets_2017_2j_1l_DeepB_eta_v2.pkl'
effFile2018ULDeepCSV = 'TTGJets_2018_2j_1l_DeepB_eta_v2.pkl'
sfFileULDeepCSV  = 'btagging.json'


effFile2016DeepJet = 'TTLep_pow_2016_2j_2l_DeepFlavB_eta_v2.pkl'
effFile2017DeepJet = 'TTLep_pow_2017_2j_2l_DeepFlavB_eta_v2.pkl'
effFile2018DeepJet = 'TTLep_pow_2018_2j_2l_DeepFlavB_eta_v2.pkl'


# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
sfFile2016DeepCSV_FastSim   = 'deepcsv_13TEV_16SL_18_3_2019.csv'
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
sfFile2017DeepCSV_FastSim   = 'deepcsv_13TEV_17SL_18_3_2019.csv'
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
sfFile2018DeepCSV_FastSim   = 'deepcsv_13TEV_18SL_7_5_2019.csv'



# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
sfFile2016DeepJet_FastSim   = 'DeepFlav_13TEV_16SL_18_3_2019.csv'
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
sfFile2017DeepJet_FastSim   = 'DeepFlav_13TEV_17SL_18_3_2019.csv'
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
sfFile2018DeepJet_FastSim   = 'DeepFlav_13TEV_18SL_7_5_2019.csv'


effFile2016preVFPULDeepJet = 'TTZToLLNuNu_UL2016_preVFP.pkl'
effFile2016postVFPULDeepJet = 'TTZToLLNuNu_UL2016.pkl'
effFile2017ULDeepJet = 'TTZToLLNuNu_UL2017.pkl'
effFile2018ULDeepJet = 'TTZToLLNuNu_UL2018.pkl'

sfFileULDeepJet  = 'btagging.json'


class BTagEfficiency:

    @staticmethod
    def getWeightDict_1b(effs, maxMultBTagWeight):
        '''Make Weight dictionary for jets
        '''
        zeroTagWeight = 1.

        for e in effs:
            zeroTagWeight*=(1-e)

        tagWeight={}
        for i in range(min(len(effs), maxMultBTagWeight)+1):
            tagWeight[i]=zeroTagWeight
            twfSum = 0.
            for tagged in itertools.combinations(effs, i):
                twf=1.
                for fac in [x/(1-x) for x in tagged]:
                    twf*=fac
                twfSum+=twf
            tagWeight[i]*=twfSum

        for i in range(maxMultBTagWeight+1):
            if not tagWeight.has_key(i):
                tagWeight[i] = 0.

        return tagWeight

    def getBTagSF_1a(self, var, bJets, nonBJets):
        if var not in self.btagWeightNames:
            raise ValueError( "Don't know what to do with b-tag variation %s" %var )
        if var != 'MC':
            ref = reduce(mul, [j['beff']['MC'] for j in bJets] + [1-j['beff']['MC'] for j in nonBJets], 1 )
            if ref>0:
                result = reduce(mul, [j['beff'][var] for j in bJets] + [1-j['beff'][var] for j in nonBJets], 1 )/ref
                return result
            else:
                logger.warning( "getBTagSF_1a: MC efficiency is zero. Return SF 1. MC efficiencies: %r "% (  [j['beff']['MC'] for j in bJets] + [1-j['beff']['MC'] for j in nonBJets] ) )
                return 1


    def __init__( self, WP=ROOT.BTagEntry.OP_MEDIUM, fastSim=False, year='UL2016', tagger='DeepCSV' ):

        if year not in [ 'UL2016', 'UL2016_preVFP', 'UL2017', 'UL2018' ]:
            raise Exception("b tag SF for year %s not known"%year)

        self.dataDir = "$CMSSW_BASE/src/Analysis/Tools/data/btagEfficiencyData/UL/"
        self.year = year
        self.tagger = tagger

        # Whether or not FS SF are to be used
        self.fastSim = fastSim

        # All btag weight names per jet
        self.btagWeightNames = [ 'MC', 'SF', 'SF_b_Down', 'SF_b_Up', 'SF_l_Down', 'SF_l_Up', 'SF_b_Down_Correlated', 'SF_b_Up_Correlated', 'SF_l_Down_Correlated', 'SF_l_Up_Correlated', 'SF_b_Down_Uncorrelated', 'SF_b_Up_Uncorrelated', 'SF_l_Down_Uncorrelated', 'SF_l_Up_Uncorrelated' ]
        if self.fastSim:
            self.btagWeightNames += [ 'SF_FS_Up', 'SF_FS_Down']

        # Input files
    	if year == 'UL2016_preVFP':
            self.etaBins = etaBins2016
    	    tag = '2016preVFP_UL'
            if tagger == 'DeepCSV':
                self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, sfFileULDeepCSV ) )
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, tag, effFile2016preVFPULDeepCSV ) )
    	    	self.WP = 'deepCSV'
            elif tagger == 'DeepJet':
                self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, sfFileULDeepJet ) )
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, tag, effFile2016preVFPULDeepJet ) )
    	    	self.WP = 'deepJet'

    	if year == 'UL2016':
            self.etaBins = etaBins2016
    	    tag = '2016postVFP_UL'
            if tagger == 'DeepCSV':
                self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, sfFileULDeepCSV ) )
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, tag, effFile2016postVFPULDeepCSV ) )
    	    	self.WP = 'deepCSV'
            elif tagger == 'DeepJet':
                self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, sfFileULDeepJet ) )
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, tag, effFile2016postVFPULDeepJet ) )
    	    	self.WP = 'deepJet'

        if year == 'UL2017':
            self.etaBins = etaBins2017
    	    tag = '2017_UL'
            if tagger == 'DeepCSV':
                self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, sfFileULDeepCSV ) )
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, tag, effFile2017ULDeepCSV ) )
    	    	self.WP = 'deepCSV'
            elif tagger == 'DeepJet':
                self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, sfFileULDeepJet ) )
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, tag, effFile2017ULDeepJet ) )
    	    	self.WP = 'deepJet'

        if year == 'UL2018':
            self.etaBins = etaBins2018
    	    tag = '2018_UL'
            if tagger == 'DeepCSV':
                self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, sfFileULDeepCSV ) )
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, tag, effFile2018ULDeepCSV ) )
    	    	self.WP = 'deepCSV'
            elif tagger == 'DeepJet':
                self.scaleFactorFile   = os.path.expandvars( os.path.join( self.dataDir, tag, sfFileULDeepJet ) )
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, tag, effFile2018ULDeepJet ) )
    	    	self.WP = 'deepJet'


        logger.info ( "Loading scale factors from %s", self.scaleFactorFile )
    	self.correction = _core.CorrectionSet.from_file(self.scaleFactorFile)

        # Load MC efficiency
        logger.info( "Loading MC efficiency %s", self.mcEfficiencyFile )
        self.mcEff = pickle.load( file( self.mcEfficiencyFile ) )

    def getMCEff(self, pdgId, pt, eta):
        ''' Get MC efficiency for jet
        '''
        for ptBin in ptBins:
            if pt>=ptBin[0] and (pt<ptBin[1] or ptBin[1]<0):
                aeta=abs(eta)
                for etaBin in self.etaBins:
                    if abs(aeta)>=etaBin[0] and abs(aeta)<etaBin[1]:
                        if abs(pdgId)==5:      return  self.mcEff[tuple(ptBin)][tuple(etaBin)]["b"]
                        elif abs(pdgId)==4:    return  self.mcEff[tuple(ptBin)][tuple(etaBin)]["c"]
                        else:                  return  self.mcEff[tuple(ptBin)][tuple(etaBin)]["other"]

        logger.debug( "No MC efficiency for pt %f eta %f pdgId %i", pt, eta, pdgId)
        return 1

    def getSF(self, pdgId, pt, eta):
        working_point = 'M'
        # BTag SF Not implemented below 20 GeV
        if pt<20:
            if self.fastSim:
                return (1,1,1,1,1,1,1)
            else:
                return (1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.)

        # BTag SF Not implemented above absEta 2.4
        if abs(eta)>=2.4:
            if self.fastSim:
                return (1,1,1,1,1,1,1)
            else:
                return (1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.)

        #autobounds are implemented now, no doubling of uncertainties necessary anymore
        flavKey = toFlavourKey(pdgId)
        #FastSim SFs
        sf_fs   = 1 #if not self.fastSim else self.readerFS.eval_auto_bounds('central', flavKey, eta, pt)
        sf_fs_u = 1 #if not self.fastSim else self.readerFS.eval_auto_bounds('down',    flavKey, eta, pt)
        sf_fs_d = 1 #if not self.fastSim else self.readerFS.eval_auto_bounds('up',      flavKey, eta, pt)
        if sf_fs == 0:  # should not happen, however, if pt=1000 (exactly) the reader will return a sf of 0.
            sf_fs = 1
            sf_fs_u = 1
            sf_fs_d = 1

        #FullSim SFs (times FSSF)
    	# UPDATE: evaluate('systematic', 'working_point', 'flavor', 'abseta', 'pt')

        if abs(pdgId)==5 or abs(pdgId)==4:
            #SF for b/c
            WP = self.WP + '_comb'
            self.evaltr = self.correction[WP]
            sf      	= sf_fs*self.evaltr.evaluate('central', working_point, abs(pdgId) , abs(eta), pt)
            sf_b_d      = sf_fs*self.evaltr.evaluate('down',    working_point, abs(pdgId) , abs(eta), pt)
            sf_b_u      = sf_fs*self.evaltr.evaluate('up',      working_point, abs(pdgId) , abs(eta), pt)
            sf_l_d  = sf
            sf_l_u  = sf
            sf_b_d_cor     = sf_fs*self.evaltr.evaluate('down_correlated',    working_point, abs(pdgId) , abs(eta), pt)
            sf_b_u_cor     = sf_fs*self.evaltr.evaluate('up_correlated',      working_point, abs(pdgId) , abs(eta), pt)
            sf_l_d_cor = sf
            sf_l_u_cor = sf
            sf_b_d_uncor     = sf_fs*self.evaltr.evaluate('down_uncorrelated',    working_point, abs(pdgId) , abs(eta), pt)
            sf_b_u_uncor     = sf_fs*self.evaltr.evaluate('up_uncorrelated',      working_point, abs(pdgId) , abs(eta), pt)
            sf_l_d_uncor = sf
            sf_l_u_uncor = sf
        else:
            #SF for light flavours
            WP = self.WP + '_incl'
            self.evaltr = self.correction[WP]
            sf      	= sf_fs*self.evaltr.evaluate('central', working_point, abs(pdgId) , abs(eta), pt)
            sf_b_d  = sf
            sf_b_u  = sf
            sf_l_d      = sf_fs*self.evaltr.evaluate('down',    working_point, abs(pdgId) , abs(eta), pt)
            sf_l_u      = sf_fs*self.evaltr.evaluate('up',      working_point, abs(pdgId) , abs(eta), pt)
            sf_b_d_cor  = sf
            sf_b_u_cor  = sf
            sf_l_d_cor      = sf_fs*self.evaltr.evaluate('down_correlated',    working_point, abs(pdgId) , abs(eta), pt)
            sf_l_u_cor      = sf_fs*self.evaltr.evaluate('up_correlated',      working_point, abs(pdgId) , abs(eta), pt)
            sf_b_d_uncor  = sf
            sf_b_u_uncor  = sf
            sf_l_d_uncor      = sf_fs*self.evaltr.evaluate('down_uncorrelated',    working_point, abs(pdgId) , abs(eta), pt)
            sf_l_u_uncor      = sf_fs*self.evaltr.evaluate('up_uncorrelated',      working_point, abs(pdgId) , abs(eta), pt)
        # print pdgId, sf, sf_b_d, sf_b_u, sf_l_d, sf_l_u
        if self.fastSim:
            return (sf, sf_b_d, sf_b_u, sf_l_d, sf_l_u, sf_b_d_cor, sf_b_u_cor, sf_l_d_cor, sf_l_u_cor, sf_b_d_uncor, sf_b_u_uncor, sf_l_d_uncor, sf_l_u_uncor, sf*sf_fs_u/sf_fs, sf*sf_fs_d/sf_fs)
        else:
            return (sf, sf_b_d, sf_b_u, sf_l_d, sf_l_u, sf_b_d_cor, sf_b_u_cor, sf_l_d_cor, sf_l_u_cor, sf_b_d_uncor, sf_b_u_uncor, sf_l_d_uncor, sf_l_u_uncor)

    def addBTagEffToJet(self, j):
        mcEff = self.getMCEff(j['hadronFlavour'], j['pt'], j['eta'])
        sf =    self.getSF(j['hadronFlavour'], j['pt'], j['eta'])
        if self.fastSim:
            j['beff'] =  {
                'MC':mcEff, 'SF':mcEff*sf[0],
                'SF_b_Down':mcEff*sf[1], 'SF_b_Up':mcEff*sf[2], 'SF_l_Down':mcEff*sf[3], 'SF_l_Up':mcEff*sf[4],
                'SF_b_Down_Correlated':mcEff*sf[5], 'SF_b_Up_Correlated':mcEff*sf[6], 'SF_l_Down_Correlated':mcEff*sf[7], 'SF_l_Up_Correlated':mcEff*sf[8],
                'SF_b_Down_Uncorrelated':mcEff*sf[9], 'SF_b_Up_Uncorrelated':mcEff*sf[10], 'SF_l_Down_Uncorrelated':mcEff*sf[11], 'SF_l_Up_Uncorrelated':mcEff*sf[12],
                'SF_FS_Up':mcEff*sf[13], 'SF_FS_Down':mcEff*sf[14],
            }
        else:
            j['beff'] =  {
                'MC':mcEff, 'SF':mcEff*sf[0],
                'SF_b_Down':mcEff*sf[1], 'SF_b_Up':mcEff*sf[2], 'SF_l_Down':mcEff*sf[3], 'SF_l_Up':mcEff*sf[4],
                'SF_b_Down_Correlated':mcEff*sf[5], 'SF_b_Up_Correlated':mcEff*sf[6], 'SF_l_Down_Correlated':mcEff*sf[7], 'SF_l_Up_Correlated':mcEff*sf[8],
                'SF_b_Down_Uncorrelated':mcEff*sf[9], 'SF_b_Up_Uncorrelated':mcEff*sf[10], 'SF_l_Down_Uncorrelated':mcEff*sf[11], 'SF_l_Up_Uncorrelated':mcEff*sf[12],
            }

#Method 1d
#https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
#https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
sfFile_1d = '$CMSSW_BASE/src/Analysis/Tools/data/btagEfficiencyData/ttH_BTV_CSVv2_13TeV_2015D_20151120.csv'
flavourSys_1d = {
    5:{'central', 'up_jes', 'down_jes', 'up_lf', 'down_lf', 'up_hfstats1', 'down_hfstats1', 'up_hfstats2', 'down_hfstats2'},
    4:{'central', 'up_cferr1', 'down_cferr1', 'up_cferr2', 'down_cferr2'},
    0:{'central', 'up_jes', 'down_jes', 'up_hf', 'down_hf', 'up_lfstats1', 'down_lfstats1', 'up_lfstats2', 'down_lfstats2'},
}
from operator import or_

class btagEfficiency_1d:

    def addBTagEffToJet(self, j):
        j['beff'] = {sys: 1. if sys not in flavourSys_1d[abs(j['hadronFlavour'])] else self.readers[sys].eval(toFlavourKey(j['hadronFlavour']), j['eta'], j['pt'], j['btagCSV']) for sys in self.btagWeightNames}

    def __init__(self,  WP = ROOT.BTagEntry.OP_MEDIUM):
        self.btagWeightNames = reduce(or_, flavourSys_1d.values())
        self.scaleFactorFile = sfFile_1d
        logger.info( "Loading scale factors from %s", self.scaleFactorFile )
        self.calib = ROOT.BTagCalibration("csvv2", self.scaleFactorFile )
        self.readers = {sys: ROOT.BTagCalibrationReader(self.calib, ROOT.BTagEntry.OP_RESHAPING, "iterativefit", sys) for sys in self.btagWeightNames}


if __name__ == "__main__":
	pass
