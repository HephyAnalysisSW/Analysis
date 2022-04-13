""" 
Extracting pre- and post-fit information, perform fits, create workspace, run fitdiagnostics, apply fit results to other data cards, and much more.


Parts stolen from:
https://github.com/HephySusySW/Workspace/blob/94X-master/DegenerateStopAnalysis/python/tools/sysTools.py
https://github.com/HephySusySW/Workspace/blob/94X-master/DegenerateStopAnalysis/python/tools/degTools.py
https://indico.cern.ch/event/747340/contributions/3198653/attachments/1744339/2823486/HComb-Tutorial-FitDiagnostics.pdf
https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/
https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/

install dependencies from combineHarvester
https://github.com/cms-analysis/CombineHarvester
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester

Initialize with the txt cardfile, but provide also the shape cardfile, the root file containing the input shapes, and if possible the workspace and fit diagnostics root file
the txt cardfile is only needed for a few additional informations, whould be nice to get rid of that
all fits and diagnostics are performed on the shape card file and workspace

naming for the input paths:
txt card file:              PATH/NAME.txt (required)
shape card file:            PATH/NAME_shapeCard.txt (not required, but should be provided)
input shapes:               PATH/NAME_shape.root (not required, but should be provided)
workspace:                  PATH/NAME_shapeCard.root (can be created with this class)
fitDiagnosticsTest:             PATH/NAME_shapeCard_FD.root (can be created with this class)
fitDiagnosticsTest (statOnly):  PATH/NAME_shapeCard_staOnly_FD.root (can be created with this class)

"""

# Standard imports
import os, sys
import math
import copy
import shutil
import uuid

import ROOT
ROOT.gROOT.SetBatch(True)

from Analysis.Tools.u_float import u_float

from RootTools.core.standard import *
from RootTools.plot.helpers  import copyIndexPHP

# Logging
import logging
logger = logging.getLogger(__name__)


class CombineResults:

    def __init__( self, cardFile, cardFileShape, plotDirectory, year, bkgOnly=False, isSearch=False, createMissingInputs=False, rebinnedCardFile=None ):
        """ cardFile: initialize with txt cardfile, but provide other cards with the proper naming
            plotDirectory: some results are copied to the plotDirectory when executing the function
            year: 2016, 2017, 2018 or 'combined'
            isSearch: False for measurements, true for searches
            createMissingInputs: if true, all needed missing data is created when initializing the class object
            rebinnedCardFile: when initializing the class object with a card from the output of CombineResults.createRebinnedResults(), provide here the card on which you performed the rebinning (fixme: get rid of that)
            e.g. newCard = CombineResults.createRebinnedResults( oldCard )
            CombineResults( cardFile = newCard, rebinnedCardFile = oldCard )
        """

        if isSearch:
            print "ATTENTION: MAY NOT WORK PROPERLY FOR SEARCHES, please check the output and provide feedback"

        if not isinstance( cardFile, str ):
            raise ValueError( "CardFile input needs to be a string with the path to the cardFile" )
        if not cardFile.endswith(".txt") or not os.path.exists( cardFile ):
            raise ValueError( "Please provide the path to the *.txt cardfile! Got: %s"%cardFile )

        self.year            = str(year)
        self.bkgOnly         = bkgOnly
        self.isSearch        = isSearch # for searches, the bkgOnly impact plots are with mu=0, for measurements mu=1
        self.txtCard         = cardFile # try to get rid of the txt cardfile!
        self.txtCardShape    = cardFileShape
        self.txtCardRebinned = rebinnedCardFile # add the original card if you created the CombineResults object with a rebinned card

        # if it is a combined card: search the txt card for information on initial cards and the directory structure
        self.channels     = self.__getChannelsFromCard( self.txtCardShape )
        

        logger.debug("channels are: {}".format(self.channels))


        # set some defaults. If a method gets some of these variables, they will be filled
        # this safes some time if they are used multiple times
        self.tRootFile           = None # needed otherwise python looses the pointer to the root file
        self.tStatOnlyRootFile   = None # needed otherwise python looses the pointer to the root file
        self.tShapeFile          = None # needed otherwise python looses the pointer to the root file
        self.binList             = None
        self.binLabels           = None
        self.processes           = None
        self.processList         = None
        self.fitResults          = None
        self.statOnlyFitResults  = None
        self.shapeInputs         = None
        self.fittedUncertainties = None
        self.constrain           = None
        self.nuisances           = None
        self.correlationHisto    = None
        self.rateParameter       = {"preFit":None, "postFit":None, "preFit_rateParam":None}
        self.estimates           = {"preFit":None, "postFit":None, "preFit_rateParam":None}
        self.uncertainties       = {"preFit":None, "postFit":None, "preFit_rateParam":None}
        self.uncertaintiesShape  = {"preFit":None, "postFit":None, "preFit_rateParam":None}
        self.pulls               = {"preFit":None, "postFit":None, "preFit_rateParam":None}
        self.covarianceHistos    = {"preFit":None, "postFit":None, "preFit_rateParam":None}
        self.regionHistos        = {"preFit":{"all":None}, "preFit_rateParam":{"all":None}, "postFit":{"all":None}}
        self.regionFile          = {"preFit":{"all":None}, "preFit_rateParam":{"all":None}, "postFit":{"all":None}}
        self.modHistos           = None



        

        #this had a different meaning in  ttgamma - try to avoid for SLS - set combinedCard as False
        logger.warning("combinedCard set to False (different to ttgamma)")
        self.combinedCard = False# len(self.channels) > 1
        self.years        = { y:y.replace("dc_","") for y in self.channels } if self.combinedCard else {self.channels[0]:self.year}
        
        # txt file for fits from shape, may be mandatory
        self.shapeCard      = cardFile.replace(".txt","_shapeCard.txt" )
        if not os.path.exists( self.shapeCard ):
            logger.warning( "Shape card not found: %s"%self.shapeCard )
            logger.warning( "Continuing with limited options!" )
            self.shapeCard = None

        # input shapes created with cardfile writer in the same step as the shapeFile, may be mandatory
        # look for information on the path to the input rootfiles in the shapeCard
        self.shapeRootFile = self.__getShapeDirectoriesFromCard( self.shapeCard ) if self.shapeCard else { dir:None for dir in self.channels }
        if not all(self.shapeRootFile.values()) or not all( map( os.path.exists, self.shapeRootFile.values() ) ):
            logger.warning( "Root file with initial shape distributions not found!" )
            logger.warning( "Continuing with limited options!" )
            self.shapeRootFile = { dir:None for dir in self.channels }

        logger.debug("shapeRootFile are: {}".format(self.shapeRootFile))
        

        # workspace from shape file
        self.rootWorkSpace = cardFile.replace(".txt","_shapeCard.root" )
        if not os.path.exists( self.rootWorkSpace ):
            if createMissingInputs:
                self.createWorkspace() # run the fit from card inputs -> result in shapeCard.root
            else:
                logger.warning( "Root card file of fit result not found: %s"%self.rootWorkSpace )
                logger.warning( "Run with create missing input to get it!" )
                self.rootWorkSpace = None


        # fit diagnostics output
        self.fitResult      = cardFile.replace(".txt","_shapeCard_FD.root" )
        if not os.path.exists( self.fitResult ):
            if createMissingInputs:
                logger.debug("creating fit diagnostics file now") 
                self.runFitDiagnostics( statOnly=False ) # run fit diagnostics
            else:
                logger.warning( "Root file of fit result not found: %s"%self.fitResult )
                logger.warning( "Continuing with limited options!" )
                self.fitResult = None

        self.fitResultStatOnly      = cardFile.replace(".txt","_shapeCard_statOnly_FD.root" )
        if not os.path.exists( self.fitResultStatOnly ):
            if createMissingInputs:
                self.runFitDiagnostics( statOnly=True ) # run stat only fit diagnostics
            else:
                logger.warning( "Root file w/ stat-only fit of fit result not found: %s"%self.fitResultStatOnly )
                logger.warning( "Continuing with limited options!" )
                self.fitResultStatOnly = None

        self.txtCardComb         = None
        self.txtCardRebinnedComb = None
        self.shapeCardComb       = None
        
        #set combinedCard artificially to False (see above)
        if self.combinedCard:
            self.txtCardComb       = {}
            self.shapeCardComb     = {}
            if rebinnedCardFile:
                self.txtCardRebinnedComb = {}
            for dir in self.channels:
                # only works for combined cards created with the cardfilewriter, as it creates sub-cards called dc_YEAR
                self.txtCardComb[dir]       = self.txtCard.replace("COMBINED", self.years[dir] )
                self.shapeCardComb[dir]     = self.shapeCard.replace("COMBINED", self.years[dir] )
                if rebinnedCardFile:
                    self.txtCardRebinnedComb[dir] = self.txtCardRebinned.replace("COMBINED", self.years[dir] )

        self.plotDirectory = plotDirectory
        if not os.path.isdir( self.plotDirectory ):
            os.makedirs( self.plotDirectory )

        

#    def __private( self ):
#    def public( self ):

    def __getChannelsFromCard( self, card ):
        """ if it is a combined card: search the txt card for information on initial cards and the directory structure
        """

        with open( card, "r" ) as f:
            imax = f.readlines()[1]

        Nchannels = int(imax.split(" ")[1]) 
        
        return ["Bin{}".format(i) for i in range(Nchannels)]


    def __getShapeDirectoriesFromCard( self, card ):
        """ look for information on the path to the input shape rootfiles in the shapeCard
        """

        with open( card, "r" ) as f:
            lines = f.readlines()

        shapes = {}
        for line in lines:
            if   not line.startswith("shapes") and not shapes.keys(): continue
            elif not line.startswith("shapes") and shapes.keys(): return shapes
            info            = line.split()
            if info[2] == "*":
                info[2] = "Bin0"
                info[3] = os.path.join( os.path.dirname( self.shapeCard ), info[3] )
            # shapes[info[2]] = info[3]
            for c in self.channels :
                shapes[c] = info[3]

        return shapes

    def __getStatOnlyFitObject( self, key=None ):
        """ get the statOnly fit objects
        """
        # return safed fitResult if available
        if self.statOnlyFitResults:
            if key and key in self.statOnlyFitResults.keys(): return self.statOnlyFitResults[key]
            elif not key:                                     return self.statOnlyFitResults

        if not self.fitResultStatOnly:
            print "Stat-only fit result not availabe, creating it!"
            self.runFitDiagnostics( statOnly=True ) # run stat only fit diagnostics

        if not self.tStatOnlyRootFile:
            self.tStatOnlyRootFile = ROOT.TFile( self.fitResultStatOnly, "READ")

        fits   = ["fit_s", "norm_prefit", "norm_fit_s", "nuisances_prefit", "nuisances_prefit_res", "shapes_prefit", "shapes_fit_s"]
        result = {}
        for fit in fits:
            result[fit] = copy.deepcopy( self.tStatOnlyRootFile.Get(fit) )
        self.statOnlyFitResults = result

        if key: return self.statOnlyFitResults[key]
        else:   return self.statOnlyFitResults

    def __getFitObject( self, key=None ):
        """ get the fit objects
        """


        print "key: {}".format(key)

        print "="*90
        print self.tRootFile
        print "="*90
        # return safed fitResult if available
        if self.fitResults:
            if key and key in self.fitResults.keys(): return self.fitResults[key]
            elif not key:                             return self.fitResults

        if not self.fitResult:
            print "Fit result not availabe, creating it!"
            self.runFitDiagnostics( statOnly=False ) # run fit diagnostics

        if not self.tRootFile:
            self.tRootFile = ROOT.TFile( self.fitResult, "READ")
            
        fits   = ["fit_b", "fit_s", "norm_prefit", "norm_fit_s", "norm_fit_b", "nuisances_prefit", "nuisances_prefit_res", "shapes_prefit", "shapes_fit_b", "shapes_fit_s"] # this is not part of our FD - yet , "overall_total_covar", "process_covar", "process_corr"]
        result = {}
        for fit in fits:
            result[fit] = copy.deepcopy( self.tRootFile.Get(fit) )
        self.fitResults = result

        if key: return self.fitResults[key]
        else:   return self.fitResults

    def __getShapeObject( self, key=None ):
        """ get the shape objects from the input shape root file
        """
        print self.shapeRootFile

        # return safed fitResult if available
        if self.shapeInputs:
            if key and key in self.shapeInputs.keys(): return self.shapeInputs[key]
            elif not key:                              return self.shapeInputs

        if not all(self.shapeRootFile.values()):
            raise ValueError( "Shape root file as input not found! Running in limited mode, thus cannot get the object needed!" )

        if not self.fitResult:
            print "Fit result not availabe, creating it!"
            self.runFitDiagnostics( statOnly=False ) # run fit diagnostics

        self.tShapeFile = {}
        result = {}

        for dir in self.channels:

            if not dir in self.tShapeFile.keys():
                self.tShapeFile[dir] = ROOT.TFile( self.shapeRootFile[dir], "READ")

            shapes = [ x.GetName() for x in self.tShapeFile[dir].GetListOfKeys() ]

            result[dir] = {}
            for shape in shapes:
                result[dir][shape] = copy.deepcopy( self.tShapeFile[dir].Get(shape) )
        self.shapeInputs = result

        if key: return self.shapeInputs[key]
        else:   return self.shapeInputs

    def __rewriteRebinnedFile( self, rootFile, postfit=False, statOnly=False, nBins=None ):
        """ rewrite the rootfile from rebinning in the style of the combine output
        """

        tRootFile = ROOT.TFile( rootFile, "READ" )

        result = {}
        if statOnly:
            fits = ["fit_s", "norm_fit_s", "nuisances_prefit", "nuisances_prefit_res"]
        else:
            fits = ["fit_b", "fit_s", "norm_fit_s", "norm_fit_b", "nuisances_prefit_res", "nuisances_prefit"]

        for fit in fits:
            result[fit] = tRootFile.Get(fit)
            try:    result[fit].SetName(fit)
            except: pass

        if statOnly:
            fits_dir   = ["shapes_prefit", "shapes_fit_s"] if postfit else ["shapes_prefit"]
        else:
            fits_dir   = ["shapes_prefit", "shapes_fit_b", "shapes_fit_s"] if postfit else ["shapes_prefit"]

        for fit in fits_dir:
            result[fit] = {}
            for d in self.channels:
                result[fit][d] = {}
                dir         = tRootFile.Get( fit+"/"+d )
                histList    = [ x.GetName() for x in dir.GetListOfKeys() if x.GetName() != "data" ] + ["data"]
                n           = nBins if nBins and nBins <= dir.Get(histList[0]).GetNbinsX() else dir.Get(histList[0]).GetNbinsX()
                # histograms have too many bins from the masked fit, remove those
                for hist in histList:
                    h = dir.Get(hist)
                    if type( h ) == ROOT.TGraphAsymmErrors:
                        dataHist = ROOT.TH1F(hist, hist, n, 0, n)
                        for i in range(n):
                            dataHist.SetBinContent(i+1, h.Eval(i+0.5))
                            dataHist.SetBinError(i+1, math.sqrt(h.Eval(i+0.5)))
                        h = dataHist.Clone()
                    elif nBins and n == nBins:
                        if hist == "total_covar":
                            covarHist = ROOT.TH2F(hist, hist, n, 0, n, n, 0, n)
                            for i in range(n):
                                for j in range(n):
                                    covarHist.SetBinContent( i+1, j+1, h.GetBinContent(i+1, j+1) )
                                    covarHist.SetBinError(   i+1, j+1, h.GetBinError(i+1, j+1)   )
                            h = covarHist.Clone()
                        else:
                            mcHist = ROOT.TH1F(hist, hist, n, 0, n)
                            for i in range(n):
                                mcHist.SetBinContent( i+1, h.GetBinContent(i+1) )
                                mcHist.SetBinError(   i+1, h.GetBinError(i+1)   )
                            h = mcHist.Clone()
                        
                    result[fit][d][hist] = copy.deepcopy(h)
            tRootFile.cd()

        tRootFile.Close()
        del tRootFile

        tRootFile = ROOT.TFile( rootFile, "RECREATE" )

        for fit in fits:
            result[fit].Write()

        for dir in fits_dir:
            tRootFile.mkdir(dir+"/")
            tRootFile.cd(dir)
            for d in self.channels:
                if not statOnly and not self.combinedCard:
                    result[dir][d]["total_covar"].SetName("process_covar")
                    result[dir][d]["total_covar"].Write()
                    result[dir][d]["total_signal"].Write()
                    result[dir][d]["total_background"].Write()
                    result[dir][d]["total"].Write()
                    result[dir][d]["total_overall"] = result[dir][d]["total"].Clone("total_overall")
                    result[dir][d]["total_overall"].SetName("total_overall")
                    result[dir][d]["total_overall"].Write()
                    result[dir][d]["data"].Write()
                tRootFile.mkdir(dir+"/"+d+"/")
                tRootFile.cd(dir+"/"+d+"/")

                for name, hist in result[dir][d].iteritems():
                    hist.Write()

                tRootFile.cd()

        tRootFile.Close()
        del tRootFile

    def __filter( self, filterDict, var=None ):
        if not isinstance( filterDict, dict ) or not var: return filterDict
        return { key:val for key, val in filterDict.items() if key==var}
            
    def __filterDict( self, filterDict, bin=None, estimate=None, nuisance=None ):
        """ filter dicts of dicts for the output
        """

        filterDict = self.__filter( filterDict, var=bin )
        if not isinstance( filterDict, dict ): return filterDict
        for b, b_dict in filterDict.items():
            filterDict[b] = self.__filter( b_dict, var=estimate )
            if not isinstance( filterDict[b], dict ): continue
            for e, e_dict in filterDict[b].items():
                filterDict[b][e] = self.__filter( e_dict, var=nuisance )
        return filterDict

    def __reduceHistogram( self, fromHisto, plotBins ):
        """ remove single bins from a histogram
        """

        newH = ROOT.TH1F( str(uuid.uuid4()), str(uuid.uuid4()), len(plotBins), 0, len(plotBins))
        j = 0
        for i in range( fromHisto.GetNbinsX() ):
            if i not in plotBins: continue
            newH.SetBinContent( j+1, fromHisto.GetBinContent(i+1) )
            newH.SetBinError(   j+1, fromHisto.GetBinError(i+1)   )
            j += 1
        self.__copyHistoSettings( fromHist=fromHisto, toHist=newH, plotBins=plotBins )
        return newH

    def __copyHistoSettings( self, fromHist, toHist, plotBins=None ):
        """ copy all our settings of a histogram when cloning it
        """

        for var in [attr for attr in dir(fromHist) if not callable(getattr(fromHist, attr)) and not attr.startswith("__")]:
            try: setattr( toHist, var, getattr( h, var ) )
            except: pass
        try:    toHist.style = fromHist.style
        except: pass
        try:    toHist.legendOption = fromHist.legendOption
        except: pass
        try:    toHist.legendText = fromHist.legendText
        except: pass

        j = 0
        for i in range( fromHist.GetNbinsX() ):
            if plotBins and i not in plotBins: continue
            toHist.GetXaxis().SetBinLabel( j+1, fromHist.GetXaxis().GetBinLabel( i+1 ) )
            j += 1

        toHist.LabelsOption("v","X")

    def __extractPOIResult( self, logFile ):
        """ extract the POI result from the logger output of the --justFit combine fit
            Simple to use (and fast) for the case where one just needs the r-value (stat-only, linearity test,...)
        """

        with open(logFile) as f:
            resultCount = 0
            for line in f:
                if 'FinalValue +/-  Error' in line: resultCount += 1
                if resultCount == 2 and '<none>' in line and line.split()[0] == "r":
                    init_val, r_stat, _, err_stat = tuple( line.split()[1:5] )
                    break
        return u_float(float(r_stat), float(err_stat))

    def __getNuisanceBinYield( self, nuisance, bin, directory, postFit=False ):
        """ get the bin yield of a +-1 sigma up variation of a certain nuisance
        """

        if directory not in self.channels:
            raise ValueError( "Directory %s unknown!"%dir )

        uncBin           = bin.replace(directory+"_","")
        yields           = self.getEstimates( postFit=postFit )[directory][uncBin]
        processes        = self.getProcessesPerBin( bin=bin )[bin]
        unc              = self.getUncertaintiesFromTxtCard(   bin=uncBin, postFit=postFit )[directory][uncBin]

        rateParamInfo    = self.getRateParameterInfo()
        rateParamInfo    = { nuisance:rateParamInfo[nuisance][bin] } if nuisance in rateParamInfo.keys() and bin in rateParamInfo[nuisance].keys() else {}
        rateParam        = { key:val.sigma/val.val if val.val else 0 for key, val in self.getRateParameter( postFit=postFit ).iteritems() }

        y, yup, ydown, sig = 0, 0, 0, 0
        # does not work for total directory, but total is not used
        for p in processes:
            # do not apply rate parameters if the process is not affected
            unc[p].update( { key:0 for key in rateParam.keys() } if not rateParamInfo or (rateParamInfo and p not in rateParamInfo[nuisance]) else rateParam )
            if p.count('signal') and self.isSearch: continue
            yproc  = yields[p].val if p in yields.keys() else 0 # yield is 0 when it is not in the results? or throw an error? FIXME
            uproc  = unc[p][nuisance]
            y     += yproc
            sig   += (yproc*uproc)**2
        yup   = y + math.sqrt(sig)
        ydown = y - math.sqrt(sig)

        return {"up":yup, "down":ydown, "relUp":yup/y if y else 0, "relDown":ydown/y if y else 0, "yield":y}


    def __regionHistos( self, postFit=False, plotBins=None, nuisances=None, bkgSubstracted=False, labelFormater=None, statOnly=False ):
        """ get a dictionary of histograms containing all region bins of the fit
            plot only certain bins with e.g. plotBins = [0, 5, 9]
            add additional +-1 sigma variation histograms of certain nuisances with 'nuisances'
            perform background substraction with bkgSubstracted=True
            each bin contains the label (from the txtcard), use the labelFormater = e.g. lambda label: label[1] to format the label in the histogram
            statOnly=True for getting the output of the statonly fit, if this fit was not performed yet, it will be then
        """

        hists    = {}
        hists_tmp = {}
        key    = "postFit" if postFit else "preFit"
        subkey = "_".join(map(str,plotBins)) if plotBins else "all"
        print subkey
        # exit(0)

        if not statOnly and subkey in self.regionHistos[key].keys() and self.regionHistos[key][subkey]:
            hists = self.regionHistos[key][subkey]
            return hists

        if   postFit and not self.bkgOnly: dirName = "shapes_fit_s"
        elif postFit and     self.bkgOnly: dirName = "shapes_fit_b"
        else:                              dirName = "shapes_prefit"
        
        print "dir name:", dirName
        print self.channels
        print "statOnly: {}".format(statOnly)
        
        if statOnly: fit = self.__getStatOnlyFitObject( key=dirName )
        else:        fit = self.__getFitObject( key=dirName )
        
        print "fit: {}".format(fit)
        print "start looping over channels"
        
        for dir in self.channels:
            # print "dir: {}".format(dir)
            # histList saves the names of objects for each Bin=Channel=CombinedRegion
            
            histList = [ x.GetName() for x in fit.Get(dir).GetListOfKeys() if x.GetName() != "data" ] + [ "data" ] 
            histList = filter( lambda hist: "total_covar" not in hist and "process_" not in hist, histList )
            histList.sort()
            
           
            
            hists[dir]    = {}
            hists_tmp[dir]    = {}
            
            for hist in histList:
                
                hists[dir][hist] = fit.Get(dir+"/"+hist).Clone()

                
                # change TGraph type to TH1F type for data
                if "data" in hist:
                    dataHist = hists[dir][histList[0]].Clone()

                    dataHist.Reset()
                    dataHist.SetName("data")

                    if type( hists[dir][hist] ) == ROOT.TGraphAsymmErrors:
                        for i in range(dataHist.GetNbinsX()):
                            print "data value: {} +/- {}".format(hists[dir][hist].Eval(i+0.5),math.sqrt(hists[dir][hist].Eval(i+0.5)))
                            dataHist.SetBinContent(i+1, hists[dir][hist].Eval(i+0.5))
                            dataHist.SetBinError(i+1, math.sqrt(hists[dir][hist].Eval(i+0.5)))
                        hists[dir]["data"] = dataHist
                    else:
                        hists[dir]["data"] = hists[dir][hist]

                    if hist != "data": del hists[dir][hist]

                    # Data Histo
                    hists[dir]["data"].style        = styles.errorStyle( ROOT.kBlack )
                    hists[dir]["data"].legendText   = "data"
                    hists[dir]["data"].legendOption = "p"


                
                
            

                if self.combinedCard:
                    k = "data" if "data" in hist else hist
                    hists[dir][k].GetXaxis().SetRangeUser(0, int(hists[dir][k].GetNbinsX()/3.))
            if nuisances:
                if isinstance( nuisances, str ): nuisances = [nuisances]
                hists[dir].update( self.getNuisanceHistosFromShapeCard( postFit=postFit, plotBins=None, nuisances=nuisances, directory=dir )[dir] )
            

            
                


            
            
        #    labels = [
        #        str(k) for k in range(0,88)
	#	      ]
	    labels = 2*(["VL","L","M","H","VH","CR"]*4+["L","M","H","VH","CR"]*4)
            # print len(labels)
	    #labels = ['allCR1aX', 'eCR1aX','muCR1aX','allCR1aY', 'eCR1aY', 'muCR1aY','allCR1bX','eCR1bX','muCR1bX','allCR1bY','eCR1bY','muCR1bY','allCR1cX','eCR1cX','muCR1cX','allCR1cY','eCR1cY','muCR1cY','allCR2aX', 'eCR2aX','muCR2aX','allCR2aY', 'eCR2aY', 'muCR2aY','allCR2bX','eCR2bX','muCR2bX','allCR2bY','eCR2bY','muCR2bY','allCR2cX','eCR2cX','muCR2cX','allCR2cY','eCR2cY','muCR2cY']
            
            if labels:
                for h_key, h in hists[dir].iteritems():
                    # print "h_key: {}".format(h_key)
                    if isinstance( h, dict ):
                        # print "h up nbins: {}".format(h["up"].GetNbinsX())
                        for i in range(h["up"].GetNbinsX()):
                            h["up"].GetXaxis().SetBinLabel( i+1, labels[i] )
                            h["down"].GetXaxis().SetBinLabel( i+1, labels[i] )
                        h["up"].LabelsOption("v","X") #"vu" for 45 degree labels
                        h["down"].LabelsOption("v","X") #"vu" for 45 degree labels
                    else:
                        # print h.GetNbinsX()
                        for i in range(h.GetNbinsX()):
                            h.GetXaxis().SetBinLabel( i+1, labels[i] )
                        h.LabelsOption("v","X") #"vu" for 45 degree labels

            # remove single bins from region plots
            if plotBins:
                for i_h, (h_key,h) in enumerate(hists[dir].iteritems()):
                    if isinstance( h, dict ):
                        hists[dir][h_key]["up"]   = self.__reduceHistogram( fromHisto=h["up"],   plotBins=plotBins )
                        hists[dir][h_key]["down"] = self.__reduceHistogram( fromHisto=h["down"], plotBins=plotBins )
                    else:
                        hists[dir][h_key] = self.__reduceHistogram( fromHisto=h, plotBins=plotBins )
        
        # print hists["Bin11"]["ZInv"].GetNbinsX() 
        
        if not statOnly: self.regionHistos[key][subkey] = hists

        if bkgSubstracted:
            for dir in self.channels:
                tot = "total" if "total" in hists[dir].keys() else "total_overall"

                # remove error on total background
                for b in range(hists[dir]["total_background"].GetNbinsX()):
                    hists[dir]["total_background"].SetBinError(b+1, 0)

                # use total - bkg as signal to get the full uncertainty
                safe = ["data","signal",tot,"total_background"]
                if nuisances: safe += nuisances
                for key in hists[dir].keys():
                    if key not in safe: del hists[dir][key]
                hists[dir]["signal"] = hists[dir][tot].Clone()

                hists[dir]["data"].Add( hists[dir]["total_background"], -1 )
                hists[dir]["signal"].Add( hists[dir]["total_background"], -1 )
                if nuisances:
                    for n in nuisances:
                        hists[dir][n]["up"].Add( hists[dir]["total_background"], -1 )
                        hists[dir][n]["down"].Add( hists[dir]["total_background"], -1 )

                hists[dir]["total_background"].Scale(0)


        return hists

    def createWorkspace( self, options="" ):
        """ do the fit and create workspace if it is not done yet
        """

        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

        if not all(self.shapeRootFile.values()):
            raise ValueError( "Shape root file as input not found! Running in limited mode, thus cannot get the object needed!" )

        cmd = "cd "+uniqueDirname+";combine --saveWorkspace -M MultiDimFit %s %s"%(options, self.shapeCard)
        print "Executing command: %s"%cmd
        os.system(cmd)

        self.rootWorkSpace = self.txtCard.replace(".txt","_shapeCard.root")
        shutil.copyfile(uniqueDirname+"/higgsCombineTest.MultiDimFit.mH120.root", self.rootWorkSpace)

        shutil.rmtree(uniqueDirname)

    def runFitDiagnostics( self, options="", statOnly=False ):
        """ run fitdiagnostics if it is not done yet
        """

        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

        if not self.rootWorkSpace:
            print "Workspace not availabe, creating it!"
            self.createWorkspace() # run the fit from card inputs

        cmd  = "cd "+uniqueDirname+";combine %s -M FitDiagnostics -t -1 --expectSignal 1 --saveNormalizations --saveWithUncertainties --saveShapes --saveOverall %s %s"%(self.rootWorkSpace, options, "--profilingMode none" if statOnly else "")
        print "Executing command: %s"%cmd
        os.system(cmd)

        if statOnly:
            self.fitResultStatOnly = self.shapeCard.replace(".txt","_statOnly_FD.root" )
            shutil.copyfile(uniqueDirname+'/fitDiagnosticsTest.root', self.fitResultStatOnly)
        else:
            self.fitResult         = self.shapeCard.replace(".txt","_FD.root" )
            shutil.copyfile(uniqueDirname+'/fitDiagnosticsTest.root', self.fitResult)

        print "Created Result%s: %s"%(" (stat-only)" if statOnly else "", self.fitResultStatOnly if statOnly else self.fitResult)
        print "Result%s: r=%s"%(" (stat-only)" if statOnly else "", self.getPulls( postFit=True, statOnly=statOnly )["r"])
        shutil.rmtree( uniqueDirname )

    def runLinearityTest( self, factor ):
        """ get the r value of a fit with scaled signal yield
            signal yield is scaled by 1/factor, the r value should then be == factor
            only a good idea for fits with expected observation (sum of MC)
        """

        print
        print "Running linearity test!"
        print "Make sure you run on expected observations!"
        print "Scaling signal events by 1./%f = %f"%(factor,1./factor)
        print

        if not all(self.shapeRootFile.values()) or not self.shapeCard:
            raise ValueError( "Input shape cards not found! Running in limited mode, thus cannot get the object needed!" )

        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

        logFile = os.path.join( uniqueDirname, "log.log" )
        shutil.copyfile( self.shapeCard, uniqueDirname+"/combineCard.txt" )
        for dir in self.channels:
            shutil.copyfile( self.shapeRootFile[dir], uniqueDirname+"/"+self.shapeRootFile[dir].split("/")[-1] )

        cmd     = "cd %s;"%uniqueDirname
        cmd    += "echo 'linTest rateParam * signal %f [0,5]' >> combineCard.txt"%(1./factor)
        print "Executing command: %s"%cmd
        os.system(cmd)

        cmd     = "cd %s;"%uniqueDirname
        cmd    += "text2workspace.py combineCard.txt;"
        cmd    += "combine --robustHesse 1 -M FitDiagnostics --justFit --setParameters linTest=%s --freezeParameters linTest -v 2 combineCard.root > log.log 2>&1"%str(1./factor)
        print "Executing command: %s"%cmd
        os.system(cmd)
        r_lin = self.__extractPOIResult( logFile )

        print "Extracted lin-test (preFit r=%f) result: %s"%(factor, r_lin)
        shutil.rmtree( uniqueDirname )

        return r_lin

    def plotPOIScan( self, rMin=0, rMax=2, points=200, addLumi=None ):
        """ https://indico.cern.ch/event/747340/contributions/3198653/attachments/1744339/2823486/HComb-Tutorial-FitDiagnostics.pdf
            create plot of likelihood with r value + splitted uncertainties (stat, syst, lumi if addLumi)
            addLumi is the part of the string that identifies all luminosity nuisances (e.g. addLumi='Luminosity')
            uncertainty is only split in stat+syst if addLumi=None
        """

        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

        if not self.rootWorkSpace:
            print "Workspace not availabe, creating it!"
            self.createWorkspace() # run the fit from card inputs

        params     = [ key for key in self.getPulls().keys() if key != "r" and not ("prop" in key and "_bin0" in key) ]
        statParams = ",".join( params )

        cmd  = "cd %s;combine -M MultiDimFit --algo grid --points %i --rMin %f --rMax %f -n bestfit --saveWorkspace %s "%(uniqueDirname,points,rMin,rMax,self.rootWorkSpace)
        cmd += ";combine -M MultiDimFit --algo grid --points %i --rMin %f --rMax %f -n stat --snapshotName MultiDimFit --freezeParameters %s higgsCombinebestfit.MultiDimFit.mH120.root"%(points,rMin,rMax,statParams)
        if addLumi:
            lumiParams = ",".join( [p for p in params if addLumi in p] )
            cmd += ";combine -M MultiDimFit --algo grid --points %i --rMin %f --rMax %f -n lumi --snapshotName MultiDimFit --freezeParameters %s higgsCombinebestfit.MultiDimFit.mH120.root"%(points,rMin,rMax,lumiParams)
            cmd += ";plot1DScan.py higgsCombinebestfit.MultiDimFit.mH120.root --output scanPOI_wLumi --others  higgsCombinelumi.MultiDimFit.mH120.root:Syst+Stat:4 higgsCombinestat.MultiDimFit.mH120.root:Stat:2 --breakdown lumi,syst,stat"
            cmd += ";mv scanPOI_wLumi.* %s/"%(self.plotDirectory)
        else:
            cmd += ";plot1DScan.py higgsCombinebestfit.MultiDimFit.mH120.root --output scanPOI --others higgsCombinestat.MultiDimFit.mH120.root:StatOnly:2 --breakdown syst,stat"
            cmd += ";mv scanPOI.* %s/"%(self.plotDirectory)
        print "Executing command: %s"%cmd
        os.system(cmd)

        print "Created POI scan plot"
        shutil.rmtree( uniqueDirname )


    def getRateParameter( self, rateParameter=None, postFit=False ):
        """ return safed rate parameter if available
        """

        key = "postFit" if postFit else "preFit"
        if self.rateParameter[key]:
            if rateParameter: return self.rateParameter[key][rateParameter]
            else:             return self.rateParameter[key]

        rateParamList = self.getRateParameterInfo().keys()
        pulls         = self.getPulls( postFit=postFit )
        rateParams    = { par:pulls[par] for par in rateParamList if par in pulls.keys() }

        self.rateParameter[key] = rateParams

        if rateParameter: return self.rateParameter[key][rateParameter]
        else:             return self.rateParameter[key]


    def getRateParameterInfo( self ):
        """ returns a dictionary with all rate parameters
            the dict contains a dict of all affected cards (e.g. for combined cards)
            which again contains a list of processes that are affected
        """

        if not self.txtCard:
            raise ValueError( "Input txt cards not found! Running in limited mode, thus cannot get the object needed!" )

        estimates = self.getProcessList( unique=True )
        rateParams = {}
        with open( self.txtCard ) as f:
            for line in f:
                if "extArg" in line or "rateParam" in line:
                    param, _, bin, proc = line.split()[:4]
                    proc = [ p for p in estimates if proc.replace("*","") in p ] if "*" in proc else [proc]
                    if param in rateParams.keys():
                        if bin in rateParams[param].keys():
                            rateParams[param][bin] += proc
                            rateParams[param][bin]  = list( set( rateParams[param][bin] ) )
                        else:
                            rateParams[param][bin] = proc
                    else:
                        rateParams[param] = { bin:proc }
        return rateParams

    def getNuisanceYields( self, nuisance, postFit=False ):
        """ return a dictionary with each bin containing the yield with the +-1sigma variation of a certain nuisance
        """
        return { dir:{ b:self.__getNuisanceBinYield( nuisance=nuisance, bin=b, directory=dir, postFit=postFit ) for b in self.getBinList( unique=True, directory=dir ) } for dir in self.channels }

    def getBinList( self, unique=True, directory=None ):
        """ get either the bin names for each process according to the cardfile ( Bin0 Bin0 Bin0 ... Bin1 Bin1 ...)
            or only the unique ones (Bin0 Bin1 ...)
            ordered list of bins
        """

        if directory and directory not in self.channels:
            raise ValueError( "Directory %s unknown!"%dir )

        # return safed binList if available
        if self.binList:
            binList = [ b for b in self.binList if directory in b ] if directory and self.combinedCard else self.binList
            if unique:
                binList = list(set(binList))
                binList = sorted(binList, key=lambda b: int(b.lower().split("bin")[-1]))
            return binList

        second = False
        with open( self.txtCard ) as f:
            for line in f:
                if len(line.split())==0: continue
                if line.split()[0].lower() == "bin" and not second:
                    second = True
                elif line.split()[0].lower() == "bin" and second:
                    binList = line.split()[1:]
                    self.binList = copy.deepcopy(binList)
                    break

        binList = [ b for b in self.binList if directory in b ] if directory and self.combinedCard else self.binList

        if unique:
            binList = list(set(binList))
            binList = sorted(binList, key=lambda b: int(b.lower().split("bin")[-1]))
            return binList

        return binList

    def tableNuisanceReport( self ):
        """ create table report of nuisance parameter ranges in the plotDirectory
        """

        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating "+uniqueDirname
        os.makedirs(uniqueDirname)
        shutil.copyfile(os.path.join(os.environ['CMSSW_BASE'], 'src', 'Analysis', 'Tools', 'python', 'cardFileWriter', 'mlfitNormsToText.py'), os.path.join(uniqueDirname, 'mlfitNormsToText.py'))

        if not self.fitResult:
            print "Fit result not availabe, creating it!"
            self.runFitDiagnostics( statOnly=False ) # run fit diagnostics

        cmd  = "cd %s;python mlfitNormsToText.py %s --uncertainties > nuisancesTable.txt"%(uniqueDirname,self.fitResult)
        cmd += ";mv nuisancesTable.txt %s/"%(self.plotDirectory)
        print "Executing command: %s"%cmd
        os.system(cmd)

        shutil.rmtree( uniqueDirname )

    def htmlNuisanceReport( self ):
        """ create html report of nuisance parameter ranges in the plotDirectory
        """

        if not self.shapeCard:
            raise ValueError( "Input shape cards not found! Running in limited mode, thus cannot get the object needed!" )

        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating "+uniqueDirname
        os.makedirs(uniqueDirname)

        shutil.copyfile(os.path.join(os.environ['CMSSW_BASE'], 'src', 'Analysis', 'Tools', 'python', 'cardFileWriter', 'systematicsAnalyzer.py'), os.path.join(uniqueDirname, 'systematicsAnalyzer.py'))

        cmd = "cd %s;python systematicsAnalyzer.py -a -f html %s > nuisanceReport.html"%(uniqueDirname, self.shapeCard)
        print "Executing command: %s"%cmd
        os.system(cmd)
        cmd = "cd %s; mv nuisanceReport.html %s/"%(uniqueDirname, self.plotDirectory)
        print "Executing command: %s"%cmd
        os.system(cmd)

        # run the same thing in "brief"
        cmd = "cd %s;python systematicsAnalyzer.py -a -f brief %s > nuisanceReport.txt"%(uniqueDirname, self.shapeCard)
        print "Executing command: %s"%cmd
        os.system(cmd)
        cmd = "cd %s; mv nuisanceReport.txt %s/"%(uniqueDirname, self.plotDirectory)
        print "Executing command: %s"%cmd
        os.system(cmd)

        shutil.rmtree( uniqueDirname )

    def printCorrelations( self, nuisance, nMax=100 ):
        """ print correlations of one nuisance to nMax other nuisances
        """

        if not self.fitResult:
            print "Fit result not availabe, creating it!"
            self.runFitDiagnostics( statOnly=False ) # run fit diagnostics

        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Printing first %i correlations of %s"%(nMax, nuisance)
        os.makedirs(uniqueDirname)
        shutil.copyfile(os.path.join(os.environ['CMSSW_BASE'], 'src', 'Analysis', 'Tools', 'python', 'cardFileWriter', 'printCorrelations.py'), os.path.join(uniqueDirname, 'printCorrelations.py'))
        cmd = "cd %s;python printCorrelations.py -i %s:fit_s -p %s --max %i"%(uniqueDirname, self.fitResult, nuisance, nMax)
        print "Executing command: %s"%cmd
        os.system(cmd)
        shutil.rmtree( uniqueDirname )

    def getBinLabels( self, labelFormater=None ):
        """ return bin labels extracted from the txt card
            labelFormater applies a function to the labels to format is accordingly
        """
        # return safed labels if available
        if self.binLabels:
            if labelFormater: 
                return { dir:map( labelFormater, self.binLabels[dir] ) for dir in self.channels }
            else: 
                return self.binLabels

        if self.combinedCard and self.txtCardRebinned: 
            cardfiles = self.txtCardRebinnedComb
        elif self.combinedCard:                       
            cardfiles = self.txtCardComb
        elif self.txtCardRebinned:                     
            cardfiles = { self.channels[0]:self.txtCardRebinned }
        else:                                          cardfiles = { self.channels[0]:self.txtCard }
        

        binLabels = {}
        for dir in [self.channels[0]]:
            binLabel = []
            with open( cardfiles[dir] ) as f:
                for line in f:
                    if line.startswith("# Bin"):
                        binLabel.append(self.years[dir] + " " + line.split(": ")[1].split("\n")[0])
                    elif line.startswith("#Muted"):
                        binLabel.append(self.years[dir] + " " + line.split(": ")[2].split("\n")[0])
            
            # print binLabel
            binLabels[dir] = binLabel

        self.binLabels = copy.copy(binLabels)
        # labelFormater   = lambda x:x.split(" ")[1].split(":")[0]
        # Bin0 all CT1>=300&&CT1<400&&HT>=400&&ISRJets_pt>=100&&l1_charge==-1&&l1_eta>=-1.5&&l1_eta<1.5&&l1_pt>=30&&mt>=0&&mt<60&&nHardBJets==0&&nSoftBJets==0

        if labelFormater:
            for dir in [self.channels[0]]:
                binLabels[dir] = map( labelFormater, binLabels[dir] )
        

        for i, bbb in enumerate(binLabels[dir]) :
            bbb[0] = "Bin{}".format(i)

        print "labels were designed"
        return binLabels

    def getNuisancesList( self, addRateParameter=True ):
        """ get a list of nuisances in the cardfile
        """
        rateParams = self.getRateParameter( postFit=True ).keys()

        if self.nuisances:
            if addRateParameter: return self.nuisances
            else:                return [ n for n in self.nuisances if n not in rateParams ]

        fit       = self.__getFitObject( key="fit_b" if self.bkgOnly else "fit_s" )
        nuisance  = fit.floatParsInit()
        nuisances = []
        iter  = nuisance.createIterator()
        var   = iter.Next()

        while var:
            if var.GetName() != "r":
                nuisances.append( var.GetName() )
            var = iter.Next()

        self.nuisances = nuisances
        if addRateParameter: return self.nuisances
        else:                return [ n for n in self.nuisances if n not in rateParams ]

    def getProcessList( self, unique=True ):
        """ get either the process names for each bin according to the cardfile ( MC1 MC2 MC3 ... MC1 MC2 ...)
            or only the unique ones (MC1 MC2 ...)
            ordered list of processes over all bins
        """

        if self.processList:
            if unique:
                processList = list(set(self.processList))
                processList.sort()
                return processList
            return self.processList

        if self.processList: return self.processList

        with open( self.txtCard ) as f:
            for line in f:
                if len(line.split())==0: continue
                if line.split()[0] == "process":
                    processList      = line.split()[1:]
                    self.processList = copy.deepcopy(processList)
                    break

        if unique:
            processList = list(set(processList))
            processList.sort()
            return processList
        return self.processList

    def getProcessesPerBin( self, bin=None ):
        """ get a dictionary with a list of processes for each bin
        """

        # return safed binList if available
        if self.processes:
            if bin: return {bin:self.processes[bin]}
            return copy.copy(self.processes)

        i           = 0
        bins        = self.getBinList( unique=True )
        procDict = {}

        with open( self.txtCard ) as f:
            for line in f:
                if len(line.split())==0: continue
                if line.split()[0] == "process":
                    if not self.processList:
                        self.processList = line.split()[1:]
                    # complex syntax needed to get the right order, set() mixes it up
                    procList = []
                    for proc in line.split()[i+1:]:
                        if proc not in procList:
                            procList.append(proc)
                        else:
                            procDict.update( {bins[i]:procList} )
                            procList = [proc]
                            i += 1
                    procDict.update( {bins[i]:procList} )
                    break
        self.processes = procDict

        if bin: return {bin:self.processes[bin]}
        return copy.copy(self.processes)

    def getPulls( self, nuisance=None, postFit=False, statOnly=False ):
        """ get a dictionary with all pulls and sigmas pre/postfit
            get the statOnly fitresults if statOnly=True, perform stat only fit if not available yet
        """

        # return safed pulls if available
        key = "postFit" if postFit else "preFit"
        # if self.pulls[key] and not statOnly:
        #     if nuisance: return self.pulls[key][nuisance]
        #     else:        return self.pulls[key]

        dirName = "fit_b" if self.bkgOnly else "fit_s"
        if statOnly: fit = self.__getStatOnlyFitObject( key=dirName )
        else:        fit = self.__getFitObject( key=dirName )
        pull  = fit.floatParsFinal() if postFit else fit.floatParsInit()
        pulls = {}
        iter  = pull.createIterator()
        var   = iter.Next()

        while var:
            pulls.update( { var.GetName():u_float(var.getValV(), var.getError()) } )
            var = iter.Next()

        if not statOnly:
            self.pulls[key] = pulls

        if nuisance: return pulls[nuisance]
        else:        return pulls


    def getUncertaintiesFromShapeCard( self, bin=None, estimate=None, nuisance=None, postFit=False ):
        """ uncertainties from shape.root card
            return safed uncertainties if available
            not yet working for flat uncertainties and rate parameters -> getUncertaintiesFromTxtCard
        """

        key = "postFit" if postFit else "preFit"
        if self.uncertaintiesShape[key]:
            if bin or estimate or nuisance:
                return {dir:self.__filterDict( self.uncertaintiesShape[key][dir], bin=bin, estimate=estimate, nuisance=nuisance ) for dir in self.channels}
            else:
                return self.uncertaintiesShape[key]

        pulls  = self.getPulls( postFit=postFit )
        allEst = self.getProcessList( unique=True )
        shapes = self.__getShapeObject()
        uncertainties = {}
        withMCStat = any( ["prop" in p for p in pulls.keys()] )
        mcStatUnc = [p for p in pulls.keys() if "prop" in p]
        
        for dir in self.channels:
            cov = self.getCovarianceHisto( directory=dir, labelFormater=None, postFit=postFit, normalize=True )
            uncertainties[dir] = {}
            for shape, shapeHisto in shapes[dir].iteritems():
                if shape.endswith("Down"): continue
                unc = [unc for unc in pulls.keys() if unc in shape and unc != "r"]
                if not unc and shape not in allEst: continue
                unc = unc[0] if unc else "stat"
                est = shape.replace("_"+unc+"Up","")
                shapeH = shapeHisto.Clone()
                if unc == "stat":
                    # stat unc
                    if not "histo" in uncertainties[dir].keys(): uncertainties[dir]["histo"] = {}
                    if not est in uncertainties[dir]["histo"].keys(): uncertainties[dir]["histo"][est] = {}
                    for i_bin in range(shapeH.GetNbinsX()):
                        _bin = "Bin%i"%i_bin
                        if not _bin in uncertainties[dir].keys(): uncertainties[dir][_bin] = {}
                        if not est in uncertainties[dir][_bin].keys(): uncertainties[dir][_bin][est] = {}
                        err = shapeH.GetBinError( i_bin+1 ) / shapeH.GetBinContent( i_bin+1 ) if shapeH.GetBinContent( i_bin+1 ) and withMCStat else 0
                        if postFit and withMCStat:
                            # log normal calculation
                            err = ((1+err)**pulls["prop_bin%s_bin%i"%(dir,i_bin)].sigma)-1
                        propH = shapeH.Clone()
                        propH.Scale(0)
                        propH.SetBinContent( i_bin+1, err )
                        shapeH.SetBinContent( i_bin+1, err )
                        uncertainties[dir][_bin][est][unc] = err
                        uncertainties[dir]["histo"][est]["prop_bin%s_bin%i"%(dir,i_bin)] = propH
                    uncertainties[dir]["histo"][est][unc] = shapeH.Clone()
                else:
                    shapeH.Add(shapes[dir][est],-1)
                    shapeH.Divide(shapes[dir][est])
                    if postFit:
                        # log normal calculation
                        for i_bin in range(shapeH.GetNbinsX()):
                            shapeH.SetBinContent(i_bin+1, ((1+shapeH.GetBinContent(i_bin+1))**pulls[unc].sigma)-1)

                    for i_bin in range(shapeH.GetNbinsX()):
                        _bin = "Bin%i"%i_bin
                        if not _bin in uncertainties[dir].keys(): uncertainties[dir][_bin] = {}
                        if not est in uncertainties[dir][_bin].keys(): uncertainties[dir][_bin][est] = {}
                        uncertainties[dir][_bin][est][unc] = shapeH.GetBinContent( i_bin+1 )
                    if not "histo" in uncertainties[dir].keys(): uncertainties[dir]["histo"] = {}
                    if not est in uncertainties[dir]["histo"].keys(): uncertainties[dir]["histo"][est] = {}
                    uncertainties[dir]["histo"][est][unc] = shapeH.Clone()
        
        self.uncertaintiesShape[key] = uncertainties

        if bin or estimate or nuisance:
            return {dir:self.__filterDict( self.uncertaintiesShape[key][dir], bin=bin, estimate=estimate, nuisance=nuisance ) for dir in self.channels}
        else:
            return self.uncertaintiesShape[key]


    def getUncertaintiesFromTxtCard( self, bin=None, estimate=None, nuisance=None, postFit=False ):
        """ uncertainties from txt card
            return safed uncertainties if available
        """

        if not self.txtCard:
            raise ValueError( "Input txt card not found! Running in limited mode, thus cannot get the object needed!" )
        
        
        
        key = "postFit" if postFit else "preFit"
        if self.uncertainties[key]:
            if bin or estimate or nuisance:
                return {dir:self.__filterDict( self.uncertainties[key][dir], bin=bin, estimate=estimate, nuisance=nuisance ) for dir in self.channels}
            else:
                return self.uncertainties[key]

        allUnc        = self.getNuisancesList( addRateParameter=False )
        allEst        = self.getProcessList( unique=False )
        pulls         = self.getPulls( postFit=postFit )
        binList       = self.getBinList( unique=False )
        uncertainties = {}

        # print allUnc
        # print allEst
        # print pulls
        # print binList 

        # exit(0)

        with open( self.txtCard ) as f:
            for line in f:
                if not line.split() or line.startswith("-"): continue
                unc = line.split()[0] 
                if unc not in allUnc: continue
                for i_bin, _bin in enumerate(binList):
                    dir = [ch for ch in self.channels if _bin.startswith(ch) or ch == "Bin0"][0]
                    _bin = _bin.replace(dir+"_","") if dir != "Bin0" else _bin
                    est = allEst[i_bin]
                    if not dir in uncertainties.keys(): uncertainties[dir] = {}
                    if not _bin in uncertainties[dir].keys(): uncertainties[dir][_bin] = {}
                    if not est in uncertainties[dir][_bin].keys(): uncertainties[dir][_bin][est] = {}
                    try:
                        uncertainties[dir][_bin][est][unc] = float(line.split()[2:][i_bin])-1
                    except:
                        uncertainties[dir][_bin][est][unc] = 0
                    if postFit and uncertainties[dir][_bin][est][unc]:
                        # log normal calculation
                        uncertainties[dir][_bin][est][unc] = ((1+uncertainties[dir][_bin][est][unc])**pulls[unc].sigma)-1

        self.uncertainties[key] = uncertainties

        if bin or estimate or nuisance:
            return {dir:self.__filterDict( self.uncertainties[key][dir], bin=bin, estimate=estimate, nuisance=nuisance ) for dir in self.channels}
        else:
            return self.uncertainties[key]

    def getObservation( self, bin=None ):
        """ return data observation in a dictionary
        """

        return {dir:{ b:b_dict["data"] for b, b_dict in o.iteritems() } for dir, o in self.getEstimates( postFit=False, bin=bin, estimate="data" ).iteritems()}

    def getEstimates( self, bin=None, estimate=None, postFit=False ):
        """ return pre/postfit estimates in a dictionary
        """

        key    = "postFit" if postFit else "preFit"
        if self.estimates[key]:
            ests = self.estimates[key]
            all = { d:self.__filterDict( dic, bin=bin, estimate=estimate ) if bin else dic for d, dic in ests.iteritems() } 
            return all

        regionHistos = self.getRegionHistos( postFit=postFit, plotBins=None )
        processes    = self.getProcessesPerBin( bin=None )
        yields       = {}
        tmp          = {}
        for dir, histoDict in regionHistos.iteritems():
            tmp[dir]     = {}
            yields[dir]  = {}
            for est, h in histoDict.iteritems():
                tmp[dir][est] = {}
                for i in range(h.GetNbinsX()):
                    y = h.GetBinContent(i+1)
                    e = h.GetBinError(i+1)
                    key = "Bin%i"%i
                    if key in tmp[dir][est].keys(): tmp[dir][est][key] += u_float( y, e )
                    else:                           tmp[dir][est][key]  = u_float( y, e )
                    yields[dir][key] = {}

            # stupid restructuring to make it compatible w/ other functions
            for b in yields[dir].keys():
                for est in tmp[dir].keys():
                    yields[dir][b][est] = tmp[dir][est][b]

        self.estimates[key] = yields

        all = { d:self.__filterDict( dic, bin=bin, estimate=estimate ) if bin else dic for d, dic in yields.iteritems() } 
        return all

    def getNuisanceHistosFromShapeCard( self, postFit=False, plotBins=None, nuisances=None, directory=None ):
        """ return a histogram of the nuisance with +- 1 sigma variation
            be careful: the quadratic sum of each postfit nuisance does not give you the full uncertainty due to (anti-)corrlations, use sumNuisanceHistos
        """
        if directory and directory not in self.channels:
            raise ValueError( "Directory %s unknown!"%dir )

        allEst   = self.getProcessList( unique=True )
        rateParams = self.getRateParameter( postFit=True ).keys()
        nuisanceHistos = {}

        dirs = self.channels if not directory else [directory]

        print dirs
        exit(0)
        for dir in dirs:
            histDict = self.getUncertaintiesFromShapeCard( postFit=postFit )[dir]["histo"]
            regions = self.getRegionHistos( postFit=postFit, plotBins=plotBins )[dir]
            nuisanceHistos[dir] = {}
            for i_n, nuisance in enumerate(nuisances):

                if nuisance in rateParams or "Lumi" in nuisance: #quick fix, use getNuisanceHistos for lnN, FIXME
                    nuisanceHistos[dir][nuisance] = self.getNuisanceHistos( postFit=postFit, plotBins=plotBins, nuisances=[nuisance], directory=dir )[dir][nuisance]
                    continue

                y = regions["total"].Clone("yield")
                y.Scale(0)
                total_err = regions["total"].Clone(nuisance)
                total_err.Scale(0)
                for est in allEst:
                    # quadratically add error histograms for each process
                    yproc = regions[est].Clone()            # process yield
                    y.Add(yproc)                            # total yield
                    if nuisance not in histDict[est].keys(): continue # nuisance does not apply to process
                    err = histDict[est][nuisance].Clone()   # relative error histogram
                    err.Multiply(yproc)                     # absolute error histogram
                    err.Multiply(err)                       # quadratically added
                    total_err.Add(err)

                for i in range(total_err.GetNbinsX()):                    
                    total_err.SetBinContent(i+1, math.sqrt(total_err.GetBinContent(i+1)))

                nuisanceHistUp   = y.Clone()
                nuisanceHistDown = y.Clone()
                nuisanceHistUp.Add( total_err )
                nuisanceHistDown.Add( total_err, -1 )

                nuisanceHistos[dir][nuisance]                    = { "up":nuisanceHistUp.Clone(), "down":nuisanceHistDown.Clone() }
                nuisanceHistos[dir][nuisance]["up"].style        = styles.lineStyle( ROOT.kSpring-1+i_n, width=3 ) #change to dynamic style
                nuisanceHistos[dir][nuisance]["down"].style      = styles.lineStyle( ROOT.kOrange+7+i_n, width=3 )
                nuisanceHistos[dir][nuisance]["up"].legendText   = nuisance + " (+1#sigma)"
                nuisanceHistos[dir][nuisance]["down"].legendText = nuisance + " (-1#sigma)"

        print nuisanceHistos
        exit(0)
        return nuisanceHistos

    def getNuisanceHistos( self, postFit=False, plotBins=None, nuisances=None, directory=None ):
        """ return a histogram of the nuisance with +- 1 sigma variation
            be careful: the quadratic sum of each postfit nuisance does not give you the full uncertainty due to (anti-)corrlations, use sumNuisanceHistos
        """

        if directory and directory not in self.channels:
            raise ValueError( "Directory %s unknown!"%dir )

        nuisanceHistos = {}

        dirs = self.channels if not directory else [directory]
        for dir in dirs:
            regions = self.getRegionHistos( postFit=postFit, plotBins=plotBins )[dir]
            nuisanceHistos[dir] = {}
            for i_n, nuisance in enumerate(nuisances):
                nuisanceHistUp      = regions["total_signal"].Clone(dir+nuisance+"up")
                nuisanceHistDown    = regions["total_signal"].Clone(dir+nuisance+"down")
                nuisanceYields  = self.getNuisanceYields( nuisance, postFit=postFit )[dir]
                for i in range(nuisanceHistUp.GetNbinsX()):
                    if self.combinedCard: key = dir+"_Bin"+str(i)
                    else:                 key = "Bin"+str(i)

                    nDict = nuisanceYields[key]
                    nuisanceHistUp.SetBinContent(   i+1, nDict["up"] )
                    nuisanceHistDown.SetBinContent( i+1, nDict["down"] )
    
                nuisanceHistos[dir][nuisance]                    = { "up":nuisanceHistUp.Clone(), "down":nuisanceHistDown.Clone() }
                nuisanceHistos[dir][nuisance]["up"].style        = styles.lineStyle( ROOT.kSpring-1+i_n, width=3 ) #change to dynamic style
                nuisanceHistos[dir][nuisance]["down"].style      = styles.lineStyle( ROOT.kOrange+7+i_n, width=3 )
                nuisanceHistos[dir][nuisance]["up"].legendText   = nuisance + " (+1#sigma)"
                nuisanceHistos[dir][nuisance]["down"].legendText = nuisance + " (-1#sigma)"

        return nuisanceHistos

    def createRebinnedResults( self, rebinningCardFile ):
        """ apply the fitresults of the current workspace on another data-card (rebinningCardFile)
            be careful: the MC stat uncertainty constraints and pulls are not applied, output is a bit more conservative (no MC stat constraints, no pulls)
            also be carful: if a nuisance of the current card is not available in the rebinningCardfile it will not be applied
            function adds the new card to the current card and performs a fit with the new card masked (not considered in the fit)
            masking applies all pulls/constraints, but does not consider the regions in the fit
            rewrite the output so that combineResults can handle it
            e.g. Results = CombineResults( cardFile=initialCard, plotDirectory=pd, year=y )
            newCard = Results.createRebinnedResults( rebinningCard )
            newResult = CombineResults( cardFile=newCard, plotDirectory=pd, year=y, rebinnedCardFile=rebinningCard )
        """

        if not self.txtCard:
            raise ValueError( "Input txt card not found! Running in limited mode, thus cannot get the object needed!" )

        # create environment
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp/", ustr)
        print "Creating "+uniqueDirname
        os.makedirs(uniqueDirname)

        resPath      = self.txtCard.replace(".txt","/")
        if not os.path.isdir(resPath): os.makedirs(resPath)
        if self.combinedCard:
            for dir in self.channels:
                if not os.path.isdir( resPath.replace("COMBINED",self.years[dir]) ): os.makedirs( resPath.replace("COMBINED",self.years[dir]) )

        optionsFit  = ""
        optionsMask = ""
        optionsTxtFit  = ""
        optionsTxtMask = ""
        mask = ""
        if self.combinedCard:
            rebinningRootFileComb = {}
            rebinningShapeCardFileComb = {}
            rebinningTxtCardFileComb = {}
            for dir in self.channels:
                # only works for combined cards created with the cardfilewriter, as it creates sub-cards called dc_YEAR
                if not os.path.isdir(resPath.replace("COMBINED",self.years[dir])): os.makedirs(resPath.replace("COMBINED",self.years[dir]))
                rebinningTxtCardFileComb[dir]   = rebinningCardFile.replace("COMBINED",self.years[dir])
                rebinningShapeCardFileComb[dir] = rebinningCardFile.replace(".txt","_shapeCard.txt").replace("COMBINED",self.years[dir])
                rebinningRootFileComb[dir]      = self.__getShapeDirectoriesFromCard( rebinningShapeCardFileComb[dir] )["Bin0"].replace("COMBINED",self.years[dir])
                optionsFit     += "fit_%s=%s "%(dir,self.shapeCardComb[dir])
                optionsMask    += "%s=%s "%(dir,rebinningShapeCardFileComb[dir])
                optionsTxtFit  += "fit_%s=%s "%(dir,self.txtCardComb[dir])
                optionsTxtMask += "%s=%s "%(dir,rebinningShapeCardFileComb[dir])
                mask           += "mask_%s=1,"%dir
            mask = mask[:-1]
        else:
            optionsFit     += "fit=%s "%(self.shapeCard)
            optionsMask    += "%s=%s "%(self.channels[0],rebinningCardFile.replace(".txt","_shapeCard.txt"))
            optionsTxtFit  += "fit=%s "%(self.txtCard)
            optionsTxtMask += "%s=%s "%(self.channels[0],rebinningCardFile)
            mask           += "mask_Bin0=1"

        f            = rebinningCardFile.split("/")[-1]
        resTxtFile   = os.path.join( resPath, f )
        resShapeFile = os.path.join( resPath, f.replace(".txt","_shapeCard.txt") )
        resShapeRoot = os.path.join( resPath, f.replace(".txt","_shape.root") )
        shutil.copyfile(os.path.join(os.environ['CMSSW_BASE'], 'src', 'Analysis', 'Tools', 'python', 'cardFileWriter', 'diffNuisances.py'), os.path.join(uniqueDirname, 'diffNuisances.py'))

        # combine fit card and muted card
        print "combining cards for muted fit"
        cmd  = "cd "+uniqueDirname+";combineCards.py %s %s > combinedCard.txt; text2workspace.py combinedCard.txt --channel-masks"%(optionsFit, optionsMask)
        cmd += ";text2workspace.py combinedCard.txt --channel-masks"
        cmd += ";combineCards.py %s %s > txtCard.txt"%(optionsTxtFit, optionsTxtMask)
        cmd += ";combine combinedCard.root -M FitDiagnostics --robustHesse 1 --forceRecreateNLL --saveShapes --saveNormalizations --saveOverall --saveWithUncertainties --setParameters %s"%mask
        print "Executing command: %s"%cmd
        os.system(cmd)

        # copy cards to final location
        logger.info("Putting combined card into dir %s", resPath)
        if self.combinedCard:
            for dir in self.channels:
                shutil.copyfile( rebinningTxtCardFileComb[dir],   resTxtFile.replace(  "COMBINED",self.years[dir]) )
                shutil.copyfile( rebinningShapeCardFileComb[dir], resShapeFile.replace("COMBINED",self.years[dir]) )
                shutil.copyfile( rebinningRootFileComb[dir],      resShapeRoot.replace("COMBINED",self.years[dir]) )

        else:
            shutil.copyfile(rebinningCardFile.replace(".txt","_shape.root"),    resShapeRoot)

        shutil.copyfile(rebinningCardFile,                                  resTxtFile)
        shutil.copyfile(rebinningCardFile.replace(".txt","_shapeCard.txt"), resShapeFile)
        shutil.copyfile(uniqueDirname+"/txtCard.txt",         resTxtFile.replace(  ".txt","_fit.txt"))
        shutil.copyfile(uniqueDirname+"/combinedCard.txt",    resShapeFile.replace(".txt","_fit.txt"))
        shutil.copyfile(uniqueDirname+"/combinedCard.root",   resShapeFile.replace(".txt",".root"))
        shutil.copyfile(uniqueDirname+"/fitDiagnosticsTest.root", resShapeFile.replace(".txt","_FD.root"))

        os.remove(uniqueDirname+"/fitDiagnosticsTest.root")

        # get number of bins of the rebinning card
        rbResults = CombineResults( cardFile=rebinningCardFile, plotDirectory=self.plotDirectory, year=self.year, bkgOnly=self.bkgOnly, isSearch=self.isSearch )
        newPulls  = rbResults.getPulls().keys()
        nBins = len( rbResults.getBinList( unique=True, directory=self.channels[0] ) )
        del rbResults

        # run fit with masked (muted) card
        print "run FitDiagnostics stat only"
        cmd = "cd "+uniqueDirname+";combine combinedCard.root --profilingMode none -M FitDiagnostics --saveWithUncertainties --saveShapes --saveNormalizations --saveOverall --setParameters %s"%mask
        print "Executing command: %s"%cmd
        os.system(cmd)

        shutil.copyfile(uniqueDirname+"/fitDiagnosticsTest.root", resShapeFile.replace(".txt","_statOnly_FD.root"))

        shutil.rmtree(uniqueDirname)

        # rewrite content in a similar way to the combine fit results
        self.__rewriteRebinnedFile( resShapeFile.replace(".txt","_FD.root"), postfit=True, nBins=nBins )
        self.__rewriteRebinnedFile( resShapeFile.replace(".txt","_statOnly_FD.root"), postfit=True, nBins=nBins, statOnly=True )

        return resTxtFile

    def getRegionHistos( self, postFit=False, plotBins=None, nuisances=None, bkgSubstracted=False, labelFormater=None, addStatOnlyHistos=False ):
        """ get a dictionary of histograms containing all region bins of the fit
            plot only certain bins with e.g. plotBins = [0, 5, 9]
            add additional +-1 sigma variation histograms of certain nuisances with 'nuisances'
            perform background substraction with bkgSubstracted=True
            each bin contains the label (from the txtcard), use the labelFormater = e.g. lambda label: label[1] to format the label in the histogram
            statOnly=True for getting the output of the statonly fit, if this fit was not performed yet, it will be then
        """

        if nuisances:
            if isinstance( nuisances, str ): nuisances = [nuisances]
            allNuisances = self.getNuisancesList() + ["totalUnc","MCStat","stat"]
            if any( [ n not in allNuisances for n in nuisances ] ):
                raise ValueError( "Any of the provided nuisances %s unknown!"%", ".join(nuisances) )
        
        hists = self.__regionHistos( postFit=postFit, plotBins=plotBins, nuisances=nuisances, bkgSubstracted=bkgSubstracted, labelFormater=labelFormater, statOnly=False )

        if addStatOnlyHistos:
            # only makes sense for postFit
            hists_stat = self.__regionHistos( postFit=True, plotBins=plotBins, nuisances=None, bkgSubstracted=bkgSubstracted, labelFormater=labelFormater, statOnly=True )
            for dir, dic in hists_stat.iteritems():
                for h_key, h in dic.iteritems():
                    hists[dir][h_key+"_stat"] = h.Clone()

        return hists

    def getRegionHistoList( self, regionHistos, processes=None, noData=False, sorted=False, bkgSubstracted=False, directory=None ):
        """ get the list of histograms and the ratio list for plotting a region plot using RootTools
            e.g.
                plots, ratioHistos = Results.getRegionHistoList( ... )
                you get: plots: [ [MCHistos], [NuisanceHistos], [dataHistos] ], ratioHistos: [ (0,0), (1,0), (2,0) ]
                plotting.draw( Plot.fromHisto( plotName, plots ), ..., ratio = { "histos":ratioHistos } )
            with sorted=True each MC in each bin will be a separated histogram sorted by the yield
        """

        print regionHistos.keys()
        print "exit"
        
        nuisances    = self.getNuisancesList() + ["totalUnc","MCStat","stat"]

        if bkgSubstracted:
            histList =  [ [ [regionHistos["signal"]] ], [(0,0)] ]
            i = 0
            for n in nuisances:
                if n in regionHistos.keys():
                    histList[0] +=  [ [regionHistos[n]["up"]], [regionHistos[n]["down"]] ]
                    histList[1] +=  [ (len(histList[0])-2,0), (len(histList[0])-1,0) ]
            histList[0] +=  [ [regionHistos["data"]] ]
            histList[1] +=  [ (len(histList[0])-1,0) ]
            return tuple(histList)

        if directory and directory not in self.channels:
            raise ValueError( "Directory %s unknown!"%dir )

        for p in processes:
            print p
            if not p in regionHistos.keys():
                # some histograms are 0, still should be in the legend
                logger.info("Histogram for %s not found! Creating one and setting it to 0! Continuing..."%p)
                regionHistos[p] = regionHistos["signal"].Clone()
                self.__copyHistoSettings( fromHist=regionHistos["signal"], toHist=regionHistos[p], plotBins=None )
                regionHistos[p].Scale(0.)
                regionHistos[p].SetName(p)
                del regionHistos[p].legendText
                regionHistos[p].notInLegend = True
    

        binProcesses = self.getProcessesPerBin()
        ratioHistos  = []
        bins         = len(self.getBinLabels()[self.channels[0]])
        i_n          = 0

        if sorted:
            histoList = [[]]
            for i in range( bins ):
                proc_list = []
                if self.combinedCard: key = directory + "_Bin%i"%i
                else:                 key = "Bin%i"%i
                for p in binProcesses[key]:
		    #for Searches, to get signal out from stack
		    if "signal" in p: continue
                    if p in regionHistos.keys():
                        # set only one bin != 0
                        tmp = regionHistos[p].Clone()
                        tmp.Scale(0.)
                        tmp.SetBinContent( i+1, regionHistos[p].GetBinContent(i+1) )
                        self.__copyHistoSettings( fromHist=regionHistos[p], toHist=tmp, plotBins=None )
                    else:
                        tmp = regionHistos["WJets"].Clone()
                        #tmp = regionHistos["signal"].Clone()
                        tmp.Scale(0.)
                        logger.info( "Adding default histogram for process %s in bin %i"%(p, i) )
                    if i != 0:
                        # remove all but the first histogram bin from the legend
                        try: del tmp.legendText
                        except: pass

                    proc_list.append(tmp)

                # sort each bin
                proc_list.sort( key=lambda h: -h.Integral() )

                # sort each bin
                histoList[0] += proc_list

        else:
            histoList = [ [p_h for p, p_h in regionHistos.iteritems() if p in binProcesses["Bin0"] and "signal" not in p ] ]
            histoList[0].sort( key=lambda h: -regionHistos[p].Integral() )

        # add data histos
        if not noData:
            histoList   += [ [regionHistos["data"]] ]
            ratioHistos += [ (1,0) ]
            i_n         += 1
	    # for searches,plot signal Histo on top, for search regions, may be not plot them
        histoList   += [ [regionHistos["signal"]] ]
	    #i_n         += 1
        # add nuisance histos at last
        for n in nuisances:
            if n in regionHistos.keys() and isinstance( regionHistos[n], dict ):
                histoList   += [ [regionHistos[n]["up"]], [regionHistos[n]["down"]] ]
                ratioHistos += [ ((i_n)*2,0),((i_n)*2+1,0) ]
                i_n         += 1

        for i in range( regionHistos["WJets"].GetNbinsX() ):
        #for i in range( regionHistos["signal"].GetNbinsX() ):
            for h_list in histoList:
                for h in h_list:
                    # make that more dynamic FIXME
                    h.GetXaxis().SetBinLabel( i+1, regionHistos["WJets"].GetXaxis().GetBinLabel( i+1 ) )
                    #h.GetXaxis().SetBinLabel( i+1, regionHistos["signal"].GetXaxis().GetBinLabel( i+1 ) )
                    h.LabelsOption("v","X") #"vu" for 45 degree labels

        return histoList, ratioHistos
        
    def getImpactPlot( self, expected=False, printPNG=False, cores=1 ):
        """ create the impact plot of the current workspace, it no workspace is provided it will be created
        """

        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

        plotName = "impacts"
        if self.bkgOnly: plotName += "_bkgOnly"
        if expected:     plotName += "_expected"

        if self.bkgOnly: options = "--freezeParameters r --setParameters r=%i"%(0 if self.isSearch else 1)
        else:            options = "--rMin 0 --rMax 10"

        if not self.rootWorkSpace:
            print "Workspace not availabe, creating it!"
            self.createWorkspace() # run the fit from card inputs

        cd             = "cd %s"%uniqueDirname
        robustFit      = "combineTool.py -M Impacts -t -1 --expectSignal 1 -m 125 -d %s --doInitialFit --robustFit 1 %s"%(self.rootWorkSpace, options)
        impactFits     = "combineTool.py -M Impacts -t -1 --expectSignal 1 -m 125 -d %s --robustFit 1 --doFits --parallel %i %s"%( self.rootWorkSpace, cores, options )
        extractImpact  = "combineTool.py -M Impacts -t -1 --expectSignal 1 -m 125 -d %s -o impacts.json"%self.rootWorkSpace
        plotImpacts    = "plotImpacts.py -i impacts.json -o %s"%plotName
        cmd            = ";".join( [ cd, robustFit, impactFits, extractImpact, plotImpacts ] )

        print "Executing command: %s"%cmd
        os.system(cmd)

        shutil.copyfile( uniqueDirname+"/%s.pdf"%plotName, "%s/%s.pdf"%(self.plotDirectory,plotName) )
        if printPNG: # useful to get a visible plot in the www directory, for nothing else
            os.system("convert -trim %s/%s.pdf -density 150 -verbose -quality 100 -flatten -sharpen 0x1.0 -geometry 1600x1600 %s/%s.png"%( self.plotDirectory, plotName, self.plotDirectory, plotName) )
            copyIndexPHP( self.plotDirectory )

        logger.info("Impact plot created at %s/%s.pdf"%(self.plotDirectory, plotName) )
        shutil.rmtree( uniqueDirname )

    def getCorrelationHisto( self ):
        """ get the correlation matrix from the current fitdiagnostics output, it no FD output is provided it will be created
        """

        if not self.fitResult:
            print "Fit result not availabe, creating it!"
            self.runFitDiagnostics( statOnly=False ) # run fit diagnostics

        if self.correlationHisto:
            return self.correlationHisto

        fit      = self.__getFitObject( key="fit_b" if self.bkgOnly else "fit_s" )
        corrhist = copy.deepcopy(fit.correlationHist())

        # bit of formating
        corrhist.GetZaxis().SetRangeUser(-1,1)
        corrhist.LabelsOption("v","X")

        self.correlationHisto = corrhist

        return self.correlationHisto

    def getCovarianceHisto( self, directory, labelFormater=None, postFit=False, normalize=False ):
        """ get the cov matrix from the current fitdiagnostics output, it no FD output is provided it will be created
        """

        if not self.fitResult:
            print "Fit result not availabe, creating it!"
            self.runFitDiagnostics( statOnly=False ) # run fit diagnostics

        key = "postFit" if postFit else "preFit"
        if self.covarianceHistos[key]:
            return self.covarianceHistos[key]

        if postFit:
            dirName = "shapes_fit_b" if self.bkgOnly else "shapes_fit_s"
        else:
            dirName = "shapes_prefit"

        fit = self.__getFitObject( key=dirName )
        matrix = copy.deepcopy( fit.Get("overall_total_covar") )
        matrix.LabelsOption("v","X")

        # set labels
        labels = self.getBinLabels( labelFormater=labelFormater )[directory]
        labels = ["Bin{}".format(i) for i in range(60)]
        # print labels
        # print matrix.GetNbinsY()
        for i in range(matrix.GetNbinsY()):
            matrix.GetYaxis().SetBinLabel( i+1, labels[i] )
            matrix.GetXaxis().SetBinLabel( i+1, labels[i] )

        self.covarianceHistos[key] = matrix.Clone()

        # normalize
        if normalize:
            norm = matrix.ProjectionX("norm")
            for i in range(matrix.GetNbinsX()):
                for j in range(matrix.GetNbinsY()):
                    mval = matrix.GetBinContent( i+1, j+1 )
                    if mval == 0 :
                        matrix.SetBinContent( i+1, j+1, 0. )
                    else :
                        matrix.SetBinContent( i+1, j+1, mval / norm.GetBinContent(i+1) )

        return matrix

    def getCorrelationMatrixEntryByNuisances( self, nuisance_i, nuisance_j ):
        """ returns the entry of the correlationMatrix according to the name of the nuisances
            function will always return 0 for correlations of mc stat with other nuisances
        """

        if nuisance_i == nuisance_j: return 1

        cov = self.getCorrelationHisto()
        index_i = None
        index_j = None

        for i in range( cov.GetNbinsX() ):
            if cov.GetXaxis().GetBinLabel( i+1 ) == nuisance_i: index_i = i+1
            if cov.GetYaxis().GetBinLabel( i+1 ) == nuisance_j: index_j = i+1
            if index_i and index_j: break

        if not (index_i and index_j):
            print "Warning: Correlation matrix entry for nuisance %s or nuisance %s not found! Returning 0!"%(nuisance_i, nuisance_j)
            return 0

        return cov.GetBinContent( index_i, index_j )



    def sumNuisanceHistos( self, hists, addStatUnc=False, postFit=False ):
        """ sum histograms using correlation matrix if postFit, sum quadratically for preFit
            hists = { "total":TH1F, Nuisance1:{"up":TH1F, "down":TH1F}, Nuisance2:{"up":TH1F, "down":TH1F} } (as in the getRegionsPlot output)
        """
        if not hists: return {}

        if not postFit and addStatUnc:
            print "Summing over preFit histograms, thus setting addStatUnc to False!"
            addStatUnc=False

        if not "total" in hists.keys():
            raise ValueError( "Total Histogram needs to be provided to sum up nuisance histos!" )

        if not "total_stat" in hists.keys() and addStatUnc:
            raise ValueError( "Total Stat Histogram not found! Please add it using getRegionHistos() with addStatOnlyHistos=True" )

        allNuisances = self.getNuisancesList()

        unKnown = [n for n in hists.keys() if not n in allNuisances + ["total"] ]
        if unKnown:
            print "Histograms not known as nuisances will not be summed: %s"%", ".join(unKnown+["total"])

        summedNuisances = [ n for n in allNuisances if n in hists.keys() ]

        if addStatUnc:
            summedNuisances += ["totalStat"]
            statUp = hists["total_stat"].Clone()
            for i in range( statUp.GetNbinsX() ):
                statUp.SetBinContent( i+1, hists["total"].GetBinContent( i+1 ) + statUp.GetBinError( i+1 ) )
            hists["totalStat"] = {"up":statUp}

        for i_n, ni in enumerate(summedNuisances):
            if "up" not in hists[ni].keys():
                raise ValueError( "Nuisance histograms not given in the format hists = { 'total':TH1F, Nuisance1:{'up':TH1F, 'down':TH1F}, Nuisance2:{'up':TH1F, 'down':TH1F} } for nuisance: %s"%ni )
            hists[ni]["up"].Add( hists["total"], -1 )

        print "Summing uncertainties of the nuisances:: %s"%", ".join(summedNuisances)

        totalRelUp = hists["total"].Clone()
        totalRelUp.Scale(0)

        for i_n, ni in enumerate(summedNuisances):
            for j_n, nj in enumerate(summedNuisances):

                if not postFit and ni != nj: continue # quadratically sum up prefit uncertainties

                tmpHist_i = hists[ni]["up"].Clone()
                tmpHist_j = hists[nj]["up"].Clone()
                tmpHist_i.Multiply( tmpHist_j )
                tmpHist_i.Scale( self.getCorrelationMatrixEntryByNuisances( ni, nj ) )
                totalRelUp.Add( tmpHist_i )

        for i in range( totalRelUp.GetNbinsX() ):
            totalRelUp.SetBinContent( i+1, math.sqrt( totalRelUp.GetBinContent( i+1 ) ) )

        totalRelDown = totalRelUp.Clone()
        totalRelDown.Scale( -1 )

        totalUp = totalRelUp.Clone()
        totalUp.Add( hists["total"] )

        totalDown = totalRelDown.Clone()
        totalDown.Add( hists["total"] )

        totalUp.legendText      = "total err (+1#sigma)"
        totalDown.legendText    = "total err (+1#sigma)"
        totalRelUp.legendText   = "total err (+1#sigma)"
        totalRelDown.legendText = "total err (+1#sigma)"

        totalUp.style        = styles.lineStyle( ROOT.kSpring-1, width=3 )
        totalRelUp.style     = styles.lineStyle( ROOT.kSpring-1, width=3 )
        totalDown.style      = styles.lineStyle( ROOT.kOrange+7, width=3 )
        totalRelDown.style   = styles.lineStyle( ROOT.kOrange+7, width=3 )

        return {"up":totalUp, "down":totalDown, "relUp":totalRelUp, "relDown":totalRelDown, "yield":hists["total"]}
