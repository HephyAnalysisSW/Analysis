#!/usr/bin/env python

# use https://github.com/danbarto/nanoAOD-tools/tree/stopsDilepton in PhysicsTools/NanoAODTools/ for this producer

# standard imports
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import os
import uuid

theano_compile_dir = '/tmp/%s'%str(uuid.uuid4())
if not os.path.exists( theano_compile_dir ):
    os.makedirs( theano_compile_dir )
os.environ['THEANO_FLAGS'] = 'base_compiledir=%s'%theano_compile_dir 

# RootTools
from RootTools.core.standard import *
from Analysis.Tools.helpers  import nonEmptyFile
### nanoAOD postprocessor
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor   import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel       import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop       import Module
    
## modules for nanoAOD postprocessor
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.METSigProducer        import METSigProducer 
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2      import *

# Logger
import logging
logger = logging.getLogger(__name__)

def extractEra(sampleName):
    return sampleName[sampleName.find("Run"):sampleName.find("Run")+len('Run2000A')]

class MetSignificance:

    def __init__( self, sample, year, output_directory, runOnUL=False ):

        if year not in [ 2016, 2017, 2018 ]:
            raise Exception("MetSignificance for year %i not known"%year)

        logger.info("Preparing nanoAOD postprocessing")
        logger.info("Will put files into directory %s", output_directory)

        self.runOnUL          = runOnUL
        self.year             = year
        self.output_directory = output_directory
        self.isData           = sample.isData
        self.name             = sample.name
        self.postfix          = "_for_%s"%self.name
        self.files            = [ f for f in sample.files if nonEmptyFile(f) ]
        self.outfiles         = None

        if self.year == 2016:
            metSigParamsMC      = [1.617529475909303, 1.617529475909303, 1.4505983036866312, 1.4505983036866312, 1.411498565372343, 1.411498565372343, 1.4087559908291813, 1.4087559908291813, 1.3633674107893856, 1.3633674107893856, 0.0019861227075085516, 0.6539410816436597]
            metSigParamsData    = [1.843242937068234, 1.843242937068234, 1.64107911184195,   1.64107911184195,   1.567040591823117, 1.567040591823117, 1.5077143780804294, 1.5077143780804294, 1.614014783345394,  1.614014783345394, -0.0005986196920895609, 0.6071479349467596]
            self.JER                 = "Summer16_25nsV1_MC" if not self.isData else "Summer16_25nsV1_DATA"
            self.jetThreshold = 15

        elif self.year == 2017:
            metSigParamsMC      = [1.9648214119268503, 1.5343086462230238, 1.9167197601498538, 1.5145044341064964, 1.8069380221985405, 1.3217263662622654, 1.5506294867561126, 1.272977540964842, 1.50742322311234,   1.6542883449796797, -0.0017865650107230548,  0.6593106706741719]
            metSigParamsData    = [2.228118299837604,  1.2420725475347338, 2.227630982417529,  1.256752205787215,  2.0215250734187853, 1.1557507029911258, 1.7350536144535336, 1.1587692458345757, 1.9385081854607988, 1.8726188460472792, -2.6697894266706265e-05, 0.646984812801919]
            self.JER                 = "Fall17_V3_MC" if not self.isData else "Fall17_V3_DATA"
            self.jetThreshold = 15

        elif self.year == 2018:
            metSigParamsMC      = [1.8430848616315363, 1.8430848616315363, 1.8572853766660877, 1.8572853766660877, 1.613083160233781,  1.613083160233781,  1.3966398718198898, 1.3966398718198898, 1.4831008506492056, 1.4831008506492056, 0.0011310724285762122, 0.6929410058142578]
            metSigParamsData    = [1.6231076732985186, 1.6231076732985186, 1.615595174619551,  1.615595174619551,  1.4731794897915416, 1.4731794897915416, 1.5183631493937553, 1.5183631493937553, 2.145670387603659,  2.145670387603659, -0.0001524158603362826, 0.7510574688006575]
            self.JER                 = "Autumn18_V7b_MC" if not self.isData else "Autumn18_V7b_DATA"
            self.jetThreshold = 25


        self.era = None
        if self.isData:
            self.era = extractEra(self.name)[-1]

        logger.info("Using JERs: %s", self.JER)

        self.unclEnThreshold = 15
        self.vetoEtaRegion   = (2.65, 3.14) if self.year == 2017 and not self.runOnUL else (10,10)
        self.METCollection   = "METFixEE2017" if self.year == 2017 and not self.runOnUL else "MET"

        # set the params for MET Significance calculation
        self.metSigParams = metSigParamsMC if not self.isData else metSigParamsData

    def __call__( self, cut ):
        newFileList = []
        logger.info("Starting nanoAOD postprocessing")
        for f in self.files:
            # need a hash to avoid data loss
            file_hash = str(hash(f))

            self.modules = []
    
            JMECorrector = createJMECorrector( isMC=(not self.isData), dataYear=self.year, runPeriod=self.era, jesUncert="Total", jetType = "AK4PFchs", metBranchName=self.METCollection, applySmearing=False )
            self.modules.append( JMECorrector() )
#            self.modules.append( METSigProducer(self.JER, self.metSigParams, METCollection=self.METCollection, useRecorr=True, calcVariations=(not self.isData), jetThreshold=self.jetThreshold, vetoEtaRegion=self.vetoEtaRegion) )

            p = PostProcessor( self.output_directory, [f], cut=cut, modules=self.modules, postfix="%s_%s"%(self.postfix, file_hash))
            p.run()
            newFileList += [self.output_directory + '/' + f.split('/')[-1].replace('.root', '%s_%s.root'%(self.postfix, file_hash))]
        logger.info("Done. Replacing input files for further processing.")
        self.outfiles = newFileList

    def getNewSampleFilenames( self ):
        # return sample.files
        return self.outfiles



if __name__ == "__main__":

    year = 2018
    mc   = False

    import os, uuid
    if year == 2016:
        if mc: from Samples.nanoAOD.Summer16_private_legacy_v1 import TTSingleLep_pow               as sample
        else:  from Samples.nanoAOD.Run2016_17Jul2018_private  import DoubleMuon_Run2016C_17Jul2018 as sample
    elif year == 2017:
        if mc: from Samples.nanoAOD.Fall17_private_legacy_v1   import TTSingleLep_pow               as sample
        else:  from Samples.nanoAOD.Run2017_31Mar2018_private  import DoubleMuon_Run2017C_31Mar2018 as sample
    elif year == 2018:
        if mc: from Samples.nanoAOD.Autumn18_private_legacy_v1 import TTSingleLep_pow               as sample
        else:  from Samples.nanoAOD.Run2018_17Sep2018_private  import DoubleMuon_Run2018C_17Sep2018 as sample

    twoJetCond             = "(Sum$(Jet_pt>=29&&abs(Jet_eta)<=2.41)>=2)"
    semilepCond_ele        = "(Sum$(Electron_pt>=34&&abs(Electron_eta)<=2.11&&Electron_cutBased>=4)>=1)"
    semilepCond_mu         = "(Sum$(Muon_pt>=29&&abs(Muon_eta)<=2.41&&Muon_tightId&&Muon_pfRelIso04_all<=0.16)>=1)"
    semilepCond            = "(" + "||".join( [semilepCond_ele, semilepCond_mu] ) + ")"
    gammaCond              = "(Sum$(Photon_pt>=19&&abs(Photon_eta)<=1.5&&Photon_electronVeto&&Photon_pixelSeed==0)>=1)"

    skimConds              = [semilepCond, gammaCond, twoJetCond]
    output_directory       = os.path.join( '/tmp/%s'%os.environ['USER'], str(uuid.uuid4()) )
    sample.reduceFiles( factor = 200 )

    MetSig = MetSignificance( sample, year, output_directory )
    MetSig( "&&".join(skimConds) )
    newfiles = MetSig.getNewSampleFilenames()
