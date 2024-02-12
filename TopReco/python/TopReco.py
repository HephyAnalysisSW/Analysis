''' Top Reconstruction interface to CMSSW module
'''

# Standard imports
import ROOT
ROOT.KinematicReconstruction
import os

# efficiciency directory

if "clip" in os.getenv("HOSTNAME").lower(): # load from CLIP
    directory = "/groups/hephy/cms/robert.schoefbeck/TopReco/data"
else:
    directory = "/afs/hephy.at/data/rschoefbeck01/TopReco/data"

def massless_LorentzVector( particle_dict ):
    return ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')( particle_dict['pt'], particle_dict['eta'], particle_dict['phi'], 0. )
    

class TopReco:

    def __init__( self, era, minNumberOfBtags, preferBtags, massLoop, b_tagger):
        self.kinReco    = ROOT.KinematicReconstruction( directory, era, minNumberOfBtags, preferBtags, massLoop)
        self.tagger     = b_tagger

    def evaluate( self, leptons, jets,  met ):
       
        if len(leptons)>=2 and leptons[0]['pdgId']*leptons[1]['pdgId']<0:
            # FIXME -> just take the leading two leptons!
            leptonMinus, leptonPlus = leptons[:2]
            if leptonMinus['pdgId']<0:
                leptonMinus, leptonPlus = leptonPlus, leptonMinus
        else:
            return

        leptonMinus_vec = massless_LorentzVector( leptonMinus )
        leptonPlus_vec  = massless_LorentzVector( leptonPlus )
        met_vec         = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(met['pt'], met['phi'], 0., 0.)
        
        jets_vec = ROOT.std.vector('ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >')()
        btags    = ROOT.std.vector('double')()
        for jet in jets:
            jets_vec.push_back( massless_LorentzVector( jet ) )
            btags.push_back( jet[ self.tagger ] )

        self.kinReco.kinReco( leptonMinus_vec, leptonPlus_vec, jets_vec, btags, met_vec )
        if self.kinReco.foundSolution: #It's important to request this, otherwise can have double counting (kinReco code returns in case of too few b-tags which lead to the same solution being returned)
            return self.kinReco.getSol() 

# 2016 top reco with >=1 btag and w/o mass loop
topReco = TopReco( ROOT.Era.run2_13tev_2016_25ns, 2, 1, 0, 'btagDeepB' )
