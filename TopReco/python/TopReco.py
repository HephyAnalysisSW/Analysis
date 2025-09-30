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

def massless_LorentzVector( particle_dict, sys_postfix=''):
    return ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')( particle_dict['pt'+sys_postfix], particle_dict['eta'], particle_dict['phi'], 0. )

class TopReco:

    def __init__( self, era, minNumberOfBtags, preferBtags, massLoop, tagger):
        self.kinReco    = ROOT.KinematicReconstruction( directory, era, minNumberOfBtags, preferBtags, massLoop)
        self.tagger = tagger

    def evaluate( self, leptons, jets,  met, sys_postfix=''):
       
        if len(leptons)>=2 and leptons[0]['pdgId']*leptons[1]['pdgId']<0:
            # FIXME -> just take the leading two leptons!
            leptonMinus, leptonPlus = leptons[:2]
            if leptonMinus['pdgId']<0:
                leptonMinus, leptonPlus = leptonPlus, leptonMinus
        else:
            return

        self.leptonMinus_vec = massless_LorentzVector( leptonMinus )
        self.leptonPlus_vec  = massless_LorentzVector( leptonPlus )
        met_vec         = ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(met['pt'], met['phi'], 0., 0.)
        
        jets_vec = ROOT.std.vector('ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >')()
        btags    = ROOT.std.vector('float')()
        for jet in jets:
            jets_vec.push_back( massless_LorentzVector( jet, sys_postfix=sys_postfix) )
            btags.push_back( jet[ self.tagger ] )

        self.kinReco.kinReco( self.leptonMinus_vec, self.leptonPlus_vec, jets_vec, btags, met_vec )
        if self.kinReco.foundSolution: #It's important to request this, otherwise can have double counting (kinReco code returns in case of too few b-tags which lead to the same solution being returned)
            sol =  self.kinReco.getSol()
            sol.leptonMinus = self.leptonMinus_vec
            sol.leptonPlus  = self.leptonPlus_vec
            return sol
        else:
            del self.leptonMinus_vec
            del self.leptonPlus_vec

if __name__=="__main__":
    # 2016 top reco with >=1 btag and w/o mass loop
    topReco = TopReco( ROOT.Era.run2_13tev_2016_25ns, 2, 1, 0, 'btagDeepB')

    sol=topReco.evaluate( leptons=[{'pt':30, 'eta':0,'phi':0,'pdgId':13}, {'pt':30, 'eta':0,'phi':3,'pdgId':-13}], jets=[{'pt':30, 'eta':0.5,'phi':0,'btagDeepB':0.99}, {'pt':30, 'eta':-0.5,'phi':0.5,'btagDeepB':0.99},{'pt':30, 'eta':-0.5,'phi':0,'btagDeepB':0.1}], met={'pt':100,'phi':0.2})
