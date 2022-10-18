''' Implements the matching of quantiles. Given two binned distributions (TH1),
    it computes x' = CDF2^-1 ( CDF1 ( x ) )
'''

# Standard imports
import ROOT
import array
import uuid

# Logger
import logging
logger = logging.getLogger(__name__)

#def transpose( h ):
#    y = [h.GetBinContent( i ) for i in range( 1, 1+h.GetNbinsX() ) ]
#    ret = ROOT.TH1D( h.GetName()+'_t', h.GetTitle()+'_t', len(y), array.array('d', y+[1]) )
#    for i in range( 1, 1+ret.GetNbinsX() ):
#        ret.SetBinContent( i, h.GetBinLowEdge(i) )
#    # overflow bin to Xmax
#    ret.SetBinContent( 1+h.GetNbinsX(), h.GetXaxis().GetXmax() ) 
#    return ret 

def transpose( h ):
    y = [ min(h.GetBinContent( i ),1) for i in range( 1, 1+h.GetNbinsX() ) ]
    name = str(uuid.uuid4())
    ret = ROOT.TH1D( name, name, len(y), array.array('d', y+[1]) )
    for i in range( 1, 1+ret.GetNbinsX() ):
        ret.SetBinContent( i, h.GetBinLowEdge(i) )
    # overflow bin to Xmax
    ret.SetBinContent( 1+h.GetNbinsX(), h.GetXaxis().GetXmax() ) 
    return ret 

def get_cumulative( h ):
    scale = h.Integral()
    h_c = h.Clone( h.GetName() + '_cumulative')
    h_c.Reset() 
    #set the overflow bin to 1
    h_c.SetBinContent( 1+h.GetNbinsX(), 1 )

    for i in range( 1, 1+h.GetNbinsX() ):
        h_c.SetBinContent(i, h_c.GetBinContent( i-1 ) + h.GetBinContent(i-1)/scale )
    return h_c

def interpolate( h, x):
    n_bins = h.GetNbinsX()
    i = h.FindBin( x )
    if i==0: return 0
    if i>n_bins: return 1
    if x-h.GetBinLowEdge(i)>0.:
        return h.GetBinContent(i) + (x-h.GetBinLowEdge(i))/h.GetBinWidth(i)*(h.GetBinContent(i+1) - h.GetBinContent(i))
    else:
        return h.GetBinContent(i)

class QuantileMatcher:

    def __init__( self, h1, h2, maxAbsU = None):
        if h1.Integral()<=0 or h2.Integral()<=0:
            raise RuntimeError( "Need positive normalization!" )

        self.h1 = h1
        self.h2 = h2

        self.maxAbsU = maxAbsU

        self.h1.Scale(1./self.h1.Integral()) 
        self.h2.Scale(1./self.h2.Integral())

        self.h1_cdf = get_cumulative(self.h1)
        self.h2_cdf = get_cumulative(self.h2)
        self.h2_cdf_inv = transpose( self.h2_cdf )

    def predict( self, u ):
        if self.maxAbsU is not None and abs(u) > self.maxAbsU:
            x = self.maxAbsU if u>0 else -self.maxAbsU 
        else:
            x = u
        
        delta = interpolate( self.h2_cdf_inv,  interpolate( self.h1_cdf,  x)) - x

        return u+delta

    def prediction_histo( self, n_bins, u_low, h_high ):
        name = str(uuid.uuid4())
        h = ROOT.TH1D(name, name,  n_bins, u_low, h_high)
        for i in range(1, n_bins+1):
            h.SetBinContent( i, self.predict( h.GetBinLowEdge(i)) )
        return h 

    @staticmethod
    def __get_quantiles( histo, quantiles):
        thresholds = array.array('d', [ROOT.Double()] * len(quantiles) )
        histo.GetQuantiles( len(quantiles), thresholds, array.array('d', quantiles) )
        return thresholds 

    def get_h1_quantiles( self, quantiles = [ 0.31, 0.5, 0.68] ):
        return self.__get_quantiles( self.h1, quantiles)
    def get_h2_quantiles( self, quantiles = [ 0.31, 0.5, 0.68] ):
        return self.__get_quantiles( self.h2, quantiles)

#if __name__=="__main__":

    #import pickle
    #data=pickle.load(file("/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/recoil_v4.3/_fine/Run2018A/lepSel-njet1p-btag0-relIso0.12-looseLeptonVeto-mll20-onZ/recoil_fitResults_SF.pkl"))
    #qm = QuantileMatcher( data[(0.0, 0.6283185307179586)][(0,50)]['para']['mc']['TH1F'], data[(0.0, 0.6283185307179586)][(0,50)]['para']['data']['TH1F'])


#    import uuid
#    ROOT.gStyle.SetOptStat(0)
#    x_range = [ -200, 200]
#    n_bins = 200
##    x_range = [ 0, 100]
##    n_bins = 100
#
#    def get_histo( f ):
#        name = str(uuid.uuid4())
#        h = ROOT.TH1D(name, name, n_bins, x_range[0],x_range[1])
#        for i in range( 1, 1+ n_bins ):
#            # Approximate probability content
#            h.SetBinContent( i, 0.5*( f.Eval( h.GetBinLowEdge( i ) ) + f.Eval( h.GetBinLowEdge(i)+h.GetBinWidth(i) )) )
#            h.SetBinError( i, 0 )
#
#        return h
#
#    f_source = ROOT.TF1('target', "[3]+[0]*exp(-0.5*((x-[1])/[2])**2)", *x_range)
#    f_source.SetLineColor(ROOT.kBlue)
#    f_source.SetParameter(0, 100)
#    f_source.SetParameter(1, 50-40)
#    f_source.SetParameter(2, 40)
#    f_source.SetParameter(3, 0)
#
#    h_source = get_histo( f_source )
#    h_source.SetLineColor(ROOT.kBlue )
#    h_source.SetLineStyle(4) 
#    h_source.Scale(1./h_source.Integral())
#
#    f_target = ROOT.TF1('target', "[3]+[0]*exp(-0.5*((x-[1])/[2])**2)", *x_range)
#    f_target.SetParameter(0, 100)
#    f_target.SetParameter(1, 50)
#    f_target.SetParameter(2, 50)
#    f_target.SetParameter(3, 0)
#
#    h_target = get_histo( f_target )
#    h_target.Scale(1./h_target.Integral())
#    h_target.SetLineStyle(4) 
#    h_target.SetLineColor(ROOT.kRed) 
#
#    # Make matcher
#    qm = QuantileMatcher( h_source, h_target )
#
#    c1 = ROOT.TCanvas()
#    #plot inputs
#    
#    c1.SetLogy(1)
#    f_source.Draw() 
#    f_target.Draw('same') 
#    c1.SetLogy()
#    c1.Print('/afs/hephy.at/user/r/rschoefbeck/www/etc/prediction_inputs.png')
#
#    c1.SetLogy(0)
#    h_source.Draw()
#    h_source.GetYaxis().SetRangeUser(0,1)
#    h_target.Draw('same')
#    qm.h1_cdf.SetLineStyle(1) 
#    qm.h2_cdf.SetLineStyle(1)
#    qm.h1_cdf.SetLineColor(ROOT.kBlue) 
#    qm.h2_cdf.SetLineColor(ROOT.kRed)
#    qm.h1_cdf.Draw('same')
#    qm.h2_cdf.Draw('same')
#    c1.Print('/afs/hephy.at/user/r/rschoefbeck/www/etc/prediction_pdf_cdf.png')
#
#    # plot prediction
#    h  = qm.prediction_histo(100,*x_range)
#    c1.SetLogy(0)
#    h.Draw()
#    h.GetYaxis().SetRangeUser(*x_range)
#    h.GetXaxis().SetRangeUser(*x_range)
#    h.GetXaxis().SetTitle( "raw u")
#    h.GetYaxis().SetTitle( "predicted u")
#    l=ROOT.TLine(x_range[0],x_range[0],x_range[1],x_range[1])
#    l.Draw()
#    c1.Print('/afs/hephy.at/user/r/rschoefbeck/www/etc/prediction_histo.png')
#
#
#    # closure
#    f_test = ROOT.TF1('test', "[3]+[0]*exp(-0.5*((x-[1])/[2])**2)", *x_range)
#    f_test.SetParameter(0, 100)
#    f_test.SetParameter(1, 50-40)
#    f_test.SetParameter(2, 50)
#    f_test.SetParameter(3, 0)
#    
#    h_test = get_histo( f_test )   
#    #h_test = h_source.Clone()
#    h_test.SetLineColor( ROOT.kBlue )
#    h_test.SetLineStyle( 1 )
#
#    h_test.Scale(1./h_test.Integral())
#    h_prediction = h_test.Clone()
#    h_prediction.Reset()
#    h_prediction.SetLineColor( ROOT.kRed )
#
#    for i in range( 1, 1+h_prediction.GetNbinsX()):
#        #h_prediction.SetBinContent( h_prediction.FindBin(qm.predict(h_test.GetBinLowEdge(i))), h_test.GetBinContent( i ) )
#        h_prediction.SetBinError( h_prediction.FindBin(qm.predict(h_test.GetBinLowEdge(i))), 0 )
#
#    h_target.Draw('h') 
#    h_source.Draw('hsame') 
#    h_test.Draw('hsame') 
#    h_prediction.Draw('hsame') 
#    c1.SetLogy(0)
#    c1.Print('/afs/hephy.at/user/r/rschoefbeck/www/etc/prediction_closure.png')
#    c1.Print('/afs/hephy.at/user/r/rschoefbeck/www/etc/prediction_closure.pdf')
