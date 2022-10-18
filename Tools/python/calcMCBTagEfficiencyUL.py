# Standard Imports
import os
import ROOT
import ROOT.TMath
import ROOT.Math
import pickle
import time
import hashlib

# Analysis Imports
from Analysis.Tools.BTagEfficiencyUL import *

def getBTagMCTruthEfficiencies( c, cut="(1)", overwrite=False, btagVar='Jet_btagCSVV2', btagWP='0.8484', etaBins=[] ):
    print(c, cut)

    etaBorders = sorted( list( set( sum( etaBins, [] ) ) ) )

    mceff = {}
    commoncf = cut + "&&"

    for ptBin in ptBins:
        mceff[ tuple(ptBin) ] = {}
        for etaBin in etaBins:
            mceff[ tuple(ptBin) ][ tuple(etaBin) ] = {}
            etaCut = "abs(Jet_eta)>" + str(etaBin[0]) + "&&abs(Jet_eta)<" + str(etaBin[1])
            ptCut  = "Jet_pt>" + str(ptBin[0])

            if ptBin[1]>0:
                ptCut += "&&Jet_pt<"+str(ptBin[1])
            c.Draw(commoncf+"("+btagVar+">"+str(btagWP)+")>>hbQuark(100,-1,2)",commoncf+"abs(Jet_hadronFlavour)==5&&                     "+etaCut+"&&"+ptCut)
            c.Draw(commoncf+"("+btagVar+">"+str(btagWP)+")>>hcQuark(100,-1,2)",commoncf+"abs(Jet_hadronFlavour)==4&&                     "+etaCut+"&&"+ptCut)
            c.Draw(commoncf+"("+btagVar+">"+str(btagWP)+")>>hOther(100,-1,2)" ,commoncf+"(abs(Jet_hadronFlavour) < 4  || abs(Jet_hadronFlavour) > 5)&&  "+etaCut+"&&"+ptCut)
            hbQuark = ROOT.gDirectory.Get("hbQuark")
            hcQuark = ROOT.gDirectory.Get("hcQuark")
            hOther = ROOT.gDirectory.Get("hOther")
            mceff[tuple(ptBin)][tuple(etaBin)]["b"]     = hbQuark.GetMean()
            mceff[tuple(ptBin)][tuple(etaBin)]["c"]     = hcQuark.GetMean()
            mceff[tuple(ptBin)][tuple(etaBin)]["other"] = hOther.GetMean()

            print("Eta",etaBin,etaCut,"Pt",ptBin,ptCut,"Found b/c/other", mceff[tuple(ptBin)][tuple(etaBin)]["b"], mceff[tuple(ptBin)][tuple(etaBin)]["c"], mceff[tuple(ptBin)][tuple(etaBin)]["other"])

            del hbQuark, hcQuark, hOther

    if overwrite: pickle.dump( mceff, file(bTagEffFile, 'w') )
    return mceff

def getBTagMCTruthEfficiencies2D( c, cut="(1)", overwrite=False, btagVar='Jet_btagCSVV2', btagWP='0.8484', etaBins=[] ):

    from array import array

    etaBorders = sorted( list( set( sum( etaBins, [] ) ) ) )

    mceff = {}
    c.SetEventList(0)

    if cut and cut.replace(" ","")!= "(1)":
        print("Setting Event List with cut: %s"%cut)
        eListName = "eList_%s"%hashlib.md5("%s"%time.time()).hexdigest()
        print(eListName)
        print(cut)
        c.Draw(">>%s"%eListName,cut)
        c.SetEventList( getattr(ROOT,eListName))

    passed_hists = {}
    total_hists = {}
    ratios = {}

    btag_var = btagVar
    btag_wp  = btagWP

    jet_quality_cut = "Jet_jetId>0"
    
    flavor_cuts = {
                        'b':'abs(Jet_hadronFlavour)==5', 
                        'c':'abs(Jet_hadronFlavour)==4',      
                        'other':'(abs(Jet_hadronFlavour) < 4  || abs(Jet_hadronFlavour) > 5)', 
                   }
   
    flavors = list(flavor_cuts.keys())
 
    for flavor in flavors:
        print("flavor", flavor)
        passed_name = 'passed_%s'%flavor
        passed_hists[flavor] = ROOT.TH2D( passed_name, passed_name , len(ptBorders)-1, array('d',ptBorders), len(etaBorders)-1, array('d', etaBorders) )
        total_name = 'total_%s'%flavor
        total_hists[flavor] = ROOT.TH2D( total_name, total_name , len(ptBorders)-1, array('d',ptBorders), len(etaBorders)-1, array('d', etaBorders) )
        c.Draw("abs(Jet_eta):Jet_pt>>%s"%passed_name, ' && '.join("(%s)"%x for x in [cut,jet_quality_cut, flavor_cuts[flavor], '%s>%s'%(btag_var, btag_wp)]))
        c.Draw("abs(Jet_eta):Jet_pt>>%s"%total_name, ' && '.join("(%s)"%x for x in [cut,jet_quality_cut, flavor_cuts[flavor] ]))
        ratios[flavor] = passed_hists[flavor].Clone("ratio_%s"%flavor)
        ratios[flavor].Divide( total_hists[flavor]) 

    for ipt, ptBin in enumerate( ptBins ,1):
        print("ptBin", ipt)
        mceff[tuple(ptBin)]={}
        for jeta, etaBin in enumerate( etaBins ,1):
            mceff[tuple(ptBin)][tuple(etaBin)] = {}
            for flavor in flavors:
                mceff[tuple(ptBin)][tuple(etaBin)][flavor] = ratios[flavor].GetBinContent(ipt, jeta)
                print("error on ptBin: %s, etaBin: %s , flavor: %s is: %s"%(ptBin,etaBin,flavor,ratios[flavor].GetBinError(ipt, jeta)))
            print("Eta",etaBin,"Pt",ptBin,"Found b/c/other", mceff[tuple(ptBin)][tuple(etaBin)]["b"], mceff[tuple(ptBin)][tuple(etaBin)]["c"], mceff[tuple(ptBin)][tuple(etaBin)]["other"])
    return mceff

def writeToFile( mcEff, filename ):
    print("write to file: ", filename)
    path = os.path.expandvars( '$CMSSW_BASE/src/Analysis/Tools/data/btagEfficiencyData/' )
    pickle.dump( mcEff, file( os.path.join( path, filename + ".pkl" ), 'w' ) )


if __name__ == "__main__":

    def get_parser():
        ''' Argument parser for post-processing module.
        '''
        import argparse
        argParser = argparse.ArgumentParser(description = "Argument parser")
        argParser.add_argument('--overwrite', action='store_true', help="Overwrite existing output files")
        argParser.add_argument('--year',      action='store', default = 'UL2016',  type=str, help="Overwrite existing output files")
        return argParser

    options = get_parser().parse_args()

    preSel = "(Sum$(Electron_pt>10&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>10&&abs(Muon_eta)<2.4)) >=3 "

    from Samples.Tools.config import redirector_global, redirector
    redirector = redirector

    if options.year == 'UL2016':
        from Samples.nanoAOD.UL16_nanoAODv9 import TTZToLLNuNu as sample
        bTagWP = '0.2598'
        etaBins = etaBins2016
    elif options.year == 'UL2016_preVFP':
        from Samples.nanoAOD.UL16_nanoAODAPVv9 import TTZToLLNuNu as sample
        bTagWP = '0.2489'
        etaBins = etaBins2016
    elif options.year == 'UL2017':
        from Samples.nanoAOD.UL17_nanoAODv9 import TTZToLLNuNu as sample
        bTagWP = '0.3040'
        etaBins = etaBins2017
    elif options.year == 'UL2018':
        from Samples.nanoAOD.UL18_nanoAODv9 import TTZToLLNuNu as sample
        bTagWP = '0.2783'
        etaBins = etaBins2018

    print(type(sample))
	
    if type(sample) == list:
        for s in sample:
            print("sample in the list", s.name)
            res = getBTagMCTruthEfficiencies2D( s.chain, cut=preSel, overwrite=options.overwrite, btagVar='Jet_btagDeepFlavB', btagWP=bTagWP, etaBins=etaBins )
            print("Efficiencies %s %s:"%(s.name, options.year))
            print(res)
            writeToFile ( res, "%s_%s"%(s.name, options.year) )
    else:
        print("sample as a chain, not list: ", sample.name)
        res = getBTagMCTruthEfficiencies2D( sample.chain, cut=preSel, overwrite=options.overwrite, btagVar='Jet_btagDeepFlavB', btagWP=bTagWP, etaBins=etaBins )
        print("Efficiencies %s %s:"%(sample.name, options.year))
        print(res)
        writeToFile ( res, "%s_%s"%(sample.name, options.year) )
    
