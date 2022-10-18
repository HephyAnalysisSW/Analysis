''' 
Extracting post-fit information.
Happily stolen from:
https://github.com/HephySusySW/Workspace/blob/94X-master/DegenerateStopAnalysis/python/tools/sysTools.py
and
https://github.com/HephySusySW/Workspace/blob/94X-master/DegenerateStopAnalysis/python/tools/degTools.py
'''

# Standard imports
import ROOT, pickle, itertools, os
import math, copy
from operator import mul

from Analysis.Tools.u_float import u_float
from Analysis.Tools.helpers import getObjFromFile

from RootTools.core.TreeReader import TreeReader

# Logging
import logging
logger = logging.getLogger(__name__)

def getValFrom1BinnedHistOrGraph( hist ):
    """
        if input is AsymTGraph, the average of errors is given 
    """
    if type(hist) in [ ROOT.TH1F , ROOT.TH1D ]:
        v = hist.GetBinContent(1)
        e = hist.GetBinError(1)
    if type(hist) in [ ROOT.TH2F , ROOT.TH2D ]:
        v = hist.GetBinContent(1,1)
        e = hist.GetBinError(1,1)
    if type(hist) in [ROOT.TGraphAsymmErrors]:
        v = hist.GetY()[0]
        el = hist.GetEYlow()[0]
        eh = hist.GetEYhigh()[0]
        if el and eh:
            e  = sum( [abs(el), abs(eh)] )/2.
        else:
            e  = max(abs(el),abs(eh))
        #print hist , (v,el,eh)
        #return (v, el, eh )
    return u_float(v,e)

def dict_function ( d,  func ):
    """
    creates a new dictionary with the same structure and depth as the input dictionary
    but the final values are determined by func(val)
    """
    new_dict = {}
    for k in list(d.keys()):
        v = d.get(k)
        if type(v)==dict:
            ret = dict_function( v , func)         
        else:
            ret = func(v)        
        new_dict[k] = ret
    return new_dict


def getPrePostFitFromMLF( mlfit ):
    if type(mlfit)==type(""):
        mlfit = ROOT.TFile(mlfit, "READ")
    shape_dirs = ['shapes_prefit', 'shapes_fit_b', 'shapes_fit_s']
    shape_hists = {}
    overalls = ['total_overall', 'total_signal', 'total_data','total_background', 'overall_total_covar'] 
    overall_outs = {}
    shape_dirs_ = {}
    for shape_dir_name in shape_dirs:
        shape_dir = mlfit.Get(shape_dir_name)
        try: shape_dir.IsFolder()
        except: continue # FIXME
        shape_dirs_[shape_dir_name]=shape_dir
        list_of_channels = [x.GetName() for x in shape_dir.GetListOfKeys() if x.IsFolder()]
        shape_hists[shape_dir_name] = {}
        overall_outs[shape_dir_name] = {}
        for channel_name in list_of_channels:
            channel  = shape_dir.Get(channel_name)
            bin_name = channel_name.replace("ch1_","")
            list_of_hists = [x.GetName() for x in channel.GetListOfKeys() ]
            shape_hists[shape_dir_name][bin_name] = {}
            for hist in list_of_hists:
                if hist =='signal' and shape_dir_name == 'shapes_fit_b' and False: ## ignore for now
                    shape_hists[shape_dir_name][bin_name][hist] = shape_dirs_['shapes_prefit'].Get(channel_name).Get(hist)
                else:
                    shape_hists[shape_dir_name][bin_name][hist] = channel.Get(hist)
                
    shape_results = dict_function( shape_hists, func = getValFrom1BinnedHistOrGraph )
    
    ret = {'hists':shape_hists, 'results':shape_results, 'mlfit':mlfit }
    return ret

def getFitResults( mlfit ):
    """ get the pre/postfit nuisances and fit results/status from the postFit root file
    """
    if type(mlfit)==type(""):
        mlfit = ROOT.TFile(mlfit, "READ")
    trees = ['tree_fit_sb', 'tree_fit_b', 'tree_prefit']
    results = {}
    for tree_name in trees:
        tree = mlfit.Get(tree_name)
        results[tree_name]={}
        list_of_channels = [x.GetName() for x in tree.GetListOfBranches()]
        for channel_name in list_of_channels:
            if "n_exp_" in channel_name or "_In" in channel_name: continue # clean up the output a bit
            tree.GetEntry(0)
            results[tree_name][channel_name] = getattr(tree, channel_name)
    return results

def getCovHisto( mlfit ):
    """ get the TH2D covariance matrix plot
    """
    if type(mlfit)==type(""):
        mlfit = ROOT.TFile(mlfit, "READ")
    shape_dirs = ['shapes_prefit', 'shapes_fit_b', 'shapes_fit_s']
    result     = {}
    for shape_dir_name in shape_dirs:
        shape_dir = mlfit.Get(shape_dir_name)
        try: shape_dir.IsFolder()
        except: continue # FIXME
        result[shape_dir_name] = copy.copy(shape_dir.Get("overall_total_covar"))

    return result

def getFitObject( mlfit ):
    """ get the fit_s and fit_b objects
    """
    if type(mlfit)==type(""):
        mlfit = ROOT.TFile(mlfit, "READ")
    fits   = ['fit_b', 'fit_s']
    result = {}
    for fit in fits:
        result[fit] = copy.copy(mlfit.Get(fit))
    return result
