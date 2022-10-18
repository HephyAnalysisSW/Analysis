''' RecoilCorrector based on QuantileMatcher
'''

# Standard imports
import ROOT
import array
import pickle

# Analysis
from Analysis.Tools.QuantileMatcher import QuantileMatcher

# Logger
import logging
logger = logging.getLogger(__name__)

class RecoilCorrector:

    def __init__(self, filename ):
        self.filename    = filename
        with open( self.filename, 'rb') as file_:
            self.correction_data = pickle.load(file_)

        self.var_bins = list(self.correction_data.keys())
        self.var_bins.sort()
        self.max_var  = max( list(map( max, self.var_bins ))) 
        self.min_var  = min( list(map( min, self.var_bins ))) 

        self.qt_bins = list(self.correction_data[self.var_bins[0]].keys())
        self.qt_bins.sort()
        self.max_qt  = max( list(map( max, self.qt_bins ))) 
        self.min_qt  = min( list(map( min, self.qt_bins ))) 

        if self.min_qt!=0:
            logger.error( "Minimum qt at %3.2f! Should be 0", self.min_qt )

        self.para_matcher = { var_bin: {qt_bin: QuantileMatcher(self.correction_data[var_bin][qt_bin]['para']['mc']['TH1F'], self.correction_data[var_bin][qt_bin]['para']['data']['TH1F']) for qt_bin in self.qt_bins} for var_bin in self.var_bins } 
        self.perp_matcher = { var_bin: {qt_bin: QuantileMatcher(self.correction_data[var_bin][qt_bin]['perp']['mc']['TH1F'], self.correction_data[var_bin][qt_bin]['perp']['data']['TH1F']) for qt_bin in self.qt_bins} for var_bin in self.var_bins } 

        logger.info( "Constructed para and perp matchers: %i var bins and %i qt bins", len(self.var_bins), len(self.qt_bins) )

    def var_bin( self, var ):
        # too low var: don't return interval
        if var<self.min_var:
            raise ValueError("Value too low.")
        # too high var: return last interval
        if var>=self.max_var:
            var = self.max_var-10**-3

        # find interval
        for iv in self.var_bins:
            if var>=iv[0] and var<iv[1]:
                return iv

    def qt_bin( self, qt ):
        # too low qt: don't return interval
        if qt<self.min_qt:
            raise ValueError("Value too low.")

        # too high qt: return last interval
        if qt>=self.max_qt:
            qt = self.max_qt-1

        # find interval
        for iv in self.qt_bins:
            if qt>=iv[0] and qt<iv[1]:
                return iv

    def predict_para(self, var, qt, u_para ):
        var_bin   = self.var_bin( var )
        qt_bin = self.qt_bin( qt )

        if var_bin and qt_bin:
            return self.para_matcher[var_bin][qt_bin].predict( u_para )

    def predict_perp(self, var, qt, u_perp ):
        var_bin   = self.var_bin( var )
        qt_bin = self.qt_bin( qt )

        if var_bin and qt_bin:
            return self.perp_matcher[var_bin][qt_bin].predict( u_perp )

if __name__=="__main__":
    import Analysis.Tools.logger as logger
    logger    = logger.get_logger( "DEBUG", logFile = None)

    #filename = '/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/recoil_v2/2018/lepSel-btag0-relIso0.12-looseLeptonVeto-mll20-onZ/recoil_fitResults_SF.pkl'
    #recoilCorrector = RecoilCorrector( filename )

