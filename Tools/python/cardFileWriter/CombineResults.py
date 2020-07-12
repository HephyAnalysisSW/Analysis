# Still to do
# signal on top of histlist for signalregions
# cardfile manipulations
# cleanup

""" 
Extracting pre- and post-fit information.
Parts stolen from:
https://github.com/HephySusySW/Workspace/blob/94X-master/DegenerateStopAnalysis/python/tools/sysTools.py
and
https://github.com/HephySusySW/Workspace/blob/94X-master/DegenerateStopAnalysis/python/tools/degTools.py

install dependencies from combineHarvester
https://github.com/cms-analysis/CombineHarvester
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester

Uses information of the *.txt, *.shapeCard.txt, shapeCard_FD.root files
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

    def __init__( self, cardFile, plotDirectory, year, bkgOnly=False, isSearch=False, createMissingInputs=False, rebinnedCardFile=None ):

        if not isinstance( cardFile, str ):
            raise ValueError( "CardFile input needs to be a string with the path to the cardFile" )
#        if not cardFile.endswith(".txt") or "shapeCard" in cardFile.split("/")[-1] or not os.path.exists( cardFile ):
        if not cardFile.endswith(".txt") or not os.path.exists( cardFile ):
            raise ValueError( "Please provide the path to the *.txt cardfile! Got: %s"%cardFile )

        self.year          = str(year)
        self.channels      = ["dc_2016","dc_2017","dc_2018"] if self.year == "combined" else ["Bin0"]
        self.bkgOnly       = bkgOnly
        self.isSearch      = isSearch # for searches, the bkgOnly impact plots are with mu=0, for measurements mu=1
        self.txtCard      = cardFile # try to get rid of the txt cardfile!
        self.txtCard16    = None
        self.txtCard17    = None
        self.txtCard18    = None
        self.txtCardRebinned = None
        self.txtCardRebinned16 = None
        self.txtCardRebinned17 = None
        self.txtCardRebinned18 = None
        if rebinnedCardFile: self.txtCardRebinned = rebinnedCardFile
        if self.year == "combined":
            self.txtCard16    = self.txtCard.replace("COMBINED","2016")
            self.txtCard17    = self.txtCard.replace("COMBINED","2017")
            self.txtCard18    = self.txtCard.replace("COMBINED","2018")
            self.shapeCard16  = None
            self.shapeCard17  = None
            self.shapeCard18  = None
            self.shapeRootFile16 = None
            self.shapeRootFile17 = None
            self.shapeRootFile18 = None
            if rebinnedCardFile:
                self.txtCardRebinned16    = self.txtCardRebinned.replace("COMBINED","2016")
                self.txtCardRebinned17    = self.txtCardRebinned.replace("COMBINED","2017")
                self.txtCardRebinned18    = self.txtCardRebinned.replace("COMBINED","2018")

        # txt file for fits from shape, may be mandatory
        self.shapeCard      = cardFile.replace(".txt","_shapeCard.txt" )
        if not os.path.exists( self.shapeCard ):
            logger.warning( "Shape card not found: %s"%self.shapeCard )
            logger.warning( "Continuing with limited options!" )
            self.shapeCard = None
        elif self.year == "combined":
            self.shapeCard16    = self.shapeCard.replace("COMBINED","2016")
            self.shapeCard17    = self.shapeCard.replace("COMBINED","2017")
            self.shapeCard18    = self.shapeCard.replace("COMBINED","2018")

        # input shapes created with cardfile writer in the same step as the shapeFile, may be mandatory
        self.shapeRootFile = cardFile.replace(".txt","_shape.root" )
        if not os.path.exists( self.shapeRootFile ):
            logger.warning( "Root file with initial shape distributions not found: %s"%self.shapeRootFile )
            logger.warning( "Continuing with limited options!" )
            self.shapeRootFile = None
        elif self.year == "combined":
            self.shapeRootFile16    = self.shapeRootFile.replace("COMBINED","2016")
            self.shapeRootFile17    = self.shapeRootFile.replace("COMBINED","2017")
            self.shapeRootFile18    = self.shapeRootFile.replace("COMBINED","2018")

        # workspace from shape file
        self.rootWorkSpace = cardFile.replace(".txt","_shapeCard.root" )
        if not os.path.exists( self.rootWorkSpace ):
            if createMissingInputs:
                self.createWorkspace() # run the fit from card inputs
            else:
                logger.warning( "Root card file of fit result not found: %s"%self.rootWorkSpace )
                logger.warning( "Creating workspace!" )
                self.rootWorkSpace = None

        # fit diagnostics output
        self.fitResult      = cardFile.replace(".txt","_shapeCard_FD.root" )
        if not os.path.exists( self.fitResult ):
            if createMissingInputs:
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

        self.plotDirectory = plotDirectory
        if not os.path.isdir( self.plotDirectory ):
            os.makedirs( self.plotDirectory )

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
        self.rateParameter       = {"preFit":None, "postFit":None}
        self.estimates           = {"preFit":None, "postFit":None}
        self.uncertainties       = {"preFit":None, "postFit":None}
        self.pulls               = {"preFit":None, "postFit":None}
        self.covarianceHistos    = {"preFit":None, "postFit":None}
        self.regionHistos        = {"preFit":{"all":None}, "postFit":{"all":None}}
        self.regionFile          = {"preFit":{"all":None}, "postFit":{"all":None}}
        self.modHistos           = None

#    def __private( self ):
#    def public( self ):

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
        # return safed fitResult if available
        if self.fitResults:
            if key and key in self.fitResults.keys(): return self.fitResults[key]
            elif not key:                             return self.fitResults

        if not self.fitResult:
            raise ValueError( "Root file of fit result not found! Running in limited mode, thus cannot get the object needed!" )

        if not self.tRootFile:
            self.tRootFile = ROOT.TFile( self.fitResult, "READ")

        fits   = ["fit_b", "fit_s", "norm_prefit", "norm_fit_s", "norm_fit_b", "nuisances_prefit", "nuisances_prefit_res", "shapes_prefit", "shapes_fit_b", "shapes_fit_s", "overall_total_covar", "process_covar", "process_corr"]
        result = {}
        for fit in fits:
            result[fit] = copy.deepcopy( self.tRootFile.Get(fit) )
        self.fitResults = result

        if key: return self.fitResults[key]
        else:   return self.fitResults

    def __getShapeObject( self, key=None ):
        """ get the fit objects
        """
        # return safed fitResult if available
        if self.shapeInputs:
            if key and key in self.shapeInputs.keys(): return self.shapeInputs[key]
            elif not key:                              return self.shapeInputs

        if not self.shapeRootFile:
            raise ValueError( "Shape root file as input not found! Running in limited mode, thus cannot get the object needed!" )

        if not self.tShapeFile:
            self.tShapeFile = ROOT.TFile( self.shapeRootFile, "READ")

        shapes = [ x.GetName() for x in self.tShapeFile.GetListOfKeys() ]

        result = {}
        for shape in shapes:
            result[shape] = copy.deepcopy( self.tShapeFile.Get(shape) )
        self.shapeInputs = result

        if key: return self.shapeInputs[key]
        else:   return self.shapeInputs

    def __rewriteRebinnedFile( self, rootFile, postfit=False, statOnly=False, nBins=None ):
        """ rewrite the rootfile from rebinning in the style of the combine output
        """
        # not yet sure what to do with combined cards
#        if self.year == "combined":
#            return

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
#            tRootFile.cd( fit )
            result[fit] = {}
            for d in self.channels:
#                tRootFile.cd( d )
                result[fit][d] = {}
                dir         = tRootFile.Get( fit+"/"+d )
                histList    = [ x.GetName() for x in dir.GetListOfKeys() if x.GetName() != "data" ] + ["data"]
                n           = nBins if nBins and nBins <= dir.Get(histList[0]).GetNbinsX() else dir.Get(histList[0]).GetNbinsX()
                # histograms have too many bins from the masked fit, remove those
                if d == "Bin0":
                    for hist in histList:
                        h = dir.Get(hist)
                        if type( h ) == ROOT.TGraphAsymmErrors:
                            dataHist = ROOT.TH1F(hist, hist, n, 0, n)
                            for i in range(n):
                                dataHist.SetBinContent(i+1, h.Eval(i+0.5))
                                dataHist.SetBinError(i+1, math.sqrt(h.Eval(i+0.5)))
                            h = dataHist.Clone()
                        elif nBins and n == nBins:
#                            h.GetXaxis().SetRangeUser(0, n)
                            if hist == "total_covar":
#                                h.GetYaxis().SetRangeUser(0, n)
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
#        tRootFile = ROOT.TFile( rootFile, "UPDATE")

        for fit in fits:
            result[fit].Write()

        for dir in fits_dir:
            tRootFile.mkdir(dir+"/")
            tRootFile.cd(dir)
            for d in self.channels:
                if not statOnly:
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
            
    def __filterDict( self, filterDict, bin=None, estimate=None, nuisance=None, systOnly=False ):
        # remove estimate sub dictionaries
        filterDict = self.__filter( filterDict, var=bin )
        if not isinstance( filterDict, dict ): return filterDict
        for b, b_dict in filterDict.items():
            filterDict[b] = self.__filter( b_dict, var=estimate )
            if not isinstance( filterDict[b], dict ): continue
            for e, e_dict in filterDict[b].items():
                filterDict[b][e] = self.__filter( e_dict, var=nuisance )
                if not isinstance( filterDict[b][e], dict ): continue
                for n, nui in filterDict[b][e].items():
                    if systOnly and ("Stat" in nui or "prop" in nui):
                        del filterDict[b][e][n]
        return filterDict

    def __reduceHistogram( self, fromHisto, plotBins ):
        # remove single bins from a histogram
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
        # copy all our settings of a histogram when cloning it
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
            # make that more dynamic FIXME
            if plotBins and i not in plotBins: continue
            toHist.GetXaxis().SetBinLabel( j+1, fromHist.GetXaxis().GetBinLabel( i+1 ) )
            j += 1

        toHist.LabelsOption("v","X")

    def __extractPOIResult( self, logFile ):
        with open(logFile) as f:
            resultCount = 0
            for line in f:
                if 'FinalValue +/-  Error' in line: resultCount += 1
                if resultCount == 2 and '<none>' in line and line.split()[0] == "r":
                    init_val, r_stat, _, err_stat = tuple( line.split()[1:5] )
                    break
        return u_float(float(r_stat), float(err_stat))#*float(init_val)

    def createWorkspace( self, options="" ):
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

        cmd = "cd "+uniqueDirname+";combine --saveWorkspace -M MultiDimFit %s %s"%(options, self.shapeCard)
        print "Executing command: %s"%cmd
        os.system(cmd)

        self.rootWorkSpace = cardFile.replace(".txt","_shapeCard.root" )
        shutil.copyfile(uniqueDirname+"/higgsCombineTest.MultiDimFit.mH120.root", self.rootWorkSpace)

        shutil.rmtree(uniqueDirname)

    def runFitDiagnostics( self, options="", statOnly=False ):
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

        cmd  = "cd "+uniqueDirname+";combine %s -M FitDiagnostics --saveNormalizations --saveWithUncertainties --saveShapes --saveOverall %s %s"%(self.rootWorkSpace, options, "--profilingMode none")
        print "Executing command: %s"%cmd
        os.system(cmd)

        if statOnly:
            self.fitResultStatOnly = self.txtCard.replace(".txt","_shapeCard_statOnly_FD.root" )
            shutil.copyfile(uniqueDirname+'/fitDiagnostics.root', self.fitResultStatOnly)
        else:
            self.fitResult         = self.txtCard.replace(".txt","_shapeCard_FD.root" )
            shutil.copyfile(uniqueDirname+'/fitDiagnostics.root', self.fitResult)

        print "Created Result%s: %s"%(" (stat-only)" if statOnly else "", self.fitResultStatOnly)
        print "Result%s: r=%s"%(" (stat-only)" if statOnly else "", self.getPulls( postFit=True, statOnly=statOnly )["r"])
        shutil.rmtree( uniqueDirname )

    def runLinearityTest( self, factor ):

        print
        print "Running linearity test!"
        print "Make sure you run on expected observations!"
        print "Scaling signal events by 1./%f = %f"%(factor,1./factor)
        print

        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

        logFile = os.path.join( uniqueDirname, "log.log" )
        shutil.copyfile( self.shapeCard, uniqueDirname+"/combineCard.txt" )
        shutil.copyfile( self.shapeRootFile, uniqueDirname+"/"+self.shapeRootFile.split("/")[-1] )

        cmd     = "cd %s;"%uniqueDirname
        cmd    += "echo 'linTest rateParam * signal %f [0,5]' >> combineCard.txt"%(1./factor)
#        cmd    += "echo 'linTest rateParam * signal' >> combineCard.txt"
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
        # https://indico.cern.ch/event/747340/contributions/3198653/attachments/1744339/2823486/HComb-Tutorial-FitDiagnostics.pdf
        # addLumi is the part of the string that identifies all luminosity nuisances (e.g. addLumi='Luminosity')
        # uncertainty is not additionally split into lumi if addLumi=None
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

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


    def getRateParameterFromTxtCard( self ):
        # returns a dictionary with all rate parameters
        # the dict contains a dict of all affected cards (e.g. for combined cards)
        # which again contains a list of processes that are affected

        binList   = self.getBinList( unique=False )
        estimates = self.getProcessList( unique=True )
        unc       = 0
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

    def getRateParameterFromShapeCard( self ):
        # not yet working
        # returns a dictionary with all rate parameters
        # the dict contains a dict of all affected cards (e.g. for combined cards)
        # which again contains a list of processes that are affected

        binList   = self.getBinList( unique=False )
        estimates = self.getProcessList( unique=True )
        unc       = 0
        rateParams = {}
        with open( self.shapeCard ) as f:
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

    def __getSubKey( self, plotBins=None ):
        # that was a bit more complicated before, now it's kinda useless
        subkey = "all"
        if plotBins: subkey = "_".join(map(str,plotBins))
        return subkey

    def getNuisanceYields( self, nuisance, postFit=False, directory="total" ):
        return { b:self.__getNuisanceBinYield( nuisance=nuisance, bin=b, postFit=postFit, directory=directory ) for b in self.getBinList( unique=True, directory=directory if self.year == "combined" and directory != "total" else None ) }

    def __getNuisanceBinYield( self, nuisance, bin, postFit=False, directory="total" ):
        yields           = self.getEstimates( postFit=postFit, directory=directory )
        processes        = self.getProcessesPerBin( bin=bin )[bin]
        unc              = self.getUncertaintiesFromTxtCard(   bin=bin, postFit=postFit, systOnly=False )[bin]

        rateParamInfo    = self.getRateParameterFromTxtCard()
        rateParamInfo    = { nuisance:rateParamInfo[nuisance][bin] } if nuisance in rateParamInfo.keys() and bin in rateParamInfo[nuisance].keys() else {}
        rateParam        = { key:val.sigma/val.val if val.val else 0 for key, val in self.getRateParameter( postFit=postFit ).iteritems() }

        if directory != "total" and self.year != "combined":
            bin = directory + "_" + bin
        elif directory == "total":
            bin = "Bin0_" + bin.split("_")[-1]

        for key, y in yields.iteritems():
            if bin in y.keys():
                yields = y[bin]
                break

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

    def getBinList( self, unique=True, directory=None ):
        # get either the bin names for each process according to the cardfile ( Bin0 Bin0 Bin0 ... Bin1 Bin1 ...)
        # or only the unique ones (Bin0 Bin1 ...)
        # ordered list of bins
        # return safed binList if available
        if self.binList:
            binList = [ b for b in self.binList if directory in b ] if directory else self.binList
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

        binList = [ b for b in self.binList if directory in b ] if directory else self.binList

        if unique:
            binList = list(set(binList))
            binList = sorted(binList, key=lambda b: int(b.lower().split("bin")[-1]))
            return binList

        return binList

    def tableNuisanceReport( self ):
        # create table report of nuisance parameter ranges
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating "+uniqueDirname
        os.makedirs(uniqueDirname)
        shutil.copyfile(os.path.join(os.environ['CMSSW_BASE'], 'src', 'Analysis', 'Tools', 'python', 'cardFileWriter', 'mlfitNormsToText.py'), os.path.join(uniqueDirname, 'mlfitNormsToText.py'))

        cmd = "cd %s;python mlfitNormsToText.py %s --uncertainties > nuisancesTable.txt"%(uniqueDirname,self.fitResult)
        cmd += ";mv nuisancesTable.txt %s/"%(self.plotDirectory)
        print "Executing command: %s"%cmd
        os.system(cmd)

        shutil.rmtree( uniqueDirname )

    def htmlNuisanceReport( self ):
        # create html report of nuisance parameter ranges
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

    def printCorrelations( self, nuisance, nMax=10 ):
        # print out correlations of nuisance
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
        # labelFormater applies a function to the labels to format is accordingly

        # return safed binList if available
        if self.binLabels:
            if labelFormater:
                binLabels = map( labelFormater, self.binLabels )
                return binLabels
            else:
                return self.binLabels

        binLabels = []
        if self.year == "combined":
            if self.txtCardRebinned:
                cardfiles = [self.txtCardRebinned16, self.txtCardRebinned17, self.txtCardRebinned18]
            else:
                cardfiles = [self.txtCard16, self.txtCard17, self.txtCard18]
            tags = ["2016 ", "2017 ", "2018 "]
        else:
            if self.txtCardRebinned:
                cardfiles = [self.txtCardRebinned]
            else:
                cardfiles = [self.txtCard]
            tags = [self.year+" "]
        for c, card in enumerate(cardfiles):
            with open( card ) as f:
                for line in f:
                    if line.startswith("# Bin"):
                        binLabels.append(tags[c] + line.split(": ")[1].split("\n")[0])
                    elif line.startswith("#Muted"):
                        binLabels.append(tags[c] + line.split(": ")[2].split("\n")[0])
        self.binLabels = binLabels

        if labelFormater:
            return map( labelFormater, binLabels )

        return self.binLabels

    def getNuisancesList( self, systOnly=False ):
        # get a list of nuisances in the cardfile

        if self.nuisances:
            if systOnly:
                for i, n in enumerate(self.nuisances):
                    if "prop" in n or "Stat" in n:
                        return self.nuisances[:i]
            return self.nuisances

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

        if systOnly:
            for i, n in enumerate(self.nuisances):
                if "prop" in n or "Stat" in n:
                    return self.nuisances[:i]

        return self.nuisances

    def getProcessList( self, unique=True ):
        # get either the process names for each bin according to the cardfile ( MC1 MC2 MC3 ... MC1 MC2 ...)
        # or only the unique ones (MC1 MC2 ...)
        # ordered list of processes over all bins
        # return safed binList if available
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
        # get a dictionary with a list of processes for each bin

        # return safed binList if available
        if self.processes:
            if bin: return {bin:self.processes[bin]}
            return self.processes

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
        return self.processes

    def getPulls( self, nuisance=None, postFit=False, statOnly=False ):

        if not self.fitResult:
            raise ValueError( "Root file of fit result not found! Running in limited mode, thus cannot get the object needed!" )

        # return safed pulls if available
        key = "postFit" if postFit else "preFit"
        if self.pulls[key] and not statOnly:
            if nuisance: return self.pulls[key][nuisance]
            else:        return self.pulls[key]

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


    def getRateParameter( self, rateParameter=None, postFit=False ):
        # return safed rate parameter if available
        key = "postFit" if postFit else "preFit"
        if self.rateParameter[key]:
            if rateParameter: return self.rateParameter[key][rateParameter]
            else:             return self.rateParameter[key]

        rateParamList = self.getRateParameterFromTxtCard().keys()
        pulls         = self.getPulls( postFit=postFit )
        rateParams    = { par:pulls[par] for par in rateParamList if par in pulls.keys() }

        self.rateParameter[key] = rateParams

        if rateParameter: return self.rateParameter[key][rateParameter]
        else:             return self.rateParameter[key]


    def getUncertaintiesFromShapeCard( self, bin=None, estimate=None, nuisance=None, postFit=False, systOnly=False ):
        # uncertainties from shape.root card
        # return safed uncertainties if available
        key = "postFit" if postFit else "preFit"
        if self.uncertainties[key]:
            if bin or estimate or nuisance or systOnly:
                return self.__filterDict( self.uncertainties[key], bin=bin, estimate=estimate, nuisance=nuisance, systOnly=systOnly )
            else:
                return self.uncertainties[key]

        pulls  = self.getPulls( postFit=postFit )
        allEst = self.getProcessList( unique=True )
        shapes = self.__getShapeObject()
        uncertainties = {}
        mcStat = any( ["prop" in p for p in pulls.keys()] )

        for shape, shapeHisto in shapes.iteritems():
            if shape.endswith("Down"): continue
            unc = [unc for unc in pulls.keys() if unc in shape and unc != "r"]
            if not unc and shape not in allEst: continue
            unc = unc[0] if unc else "stat"
            est = shape.replace("_"+unc+"Up","")
            shapeH = shapeHisto.Clone()
            if unc == "stat":
                # stat unc
                for i_bin in range(shapeH.GetNbinsX()):
                    _bin = "Bin"+str(i_bin)
                    if not _bin in uncertainties.keys(): uncertainties[_bin] = {}
                    if not est in uncertainties[_bin].keys(): uncertainties[_bin][est] = {}
                    err = shapeH.GetBinError( i_bin+1 ) / shapeH.GetBinContent( i_bin+1 ) if shapeH.GetBinContent( i_bin+1 ) and mcStat else 0
                    if postFit and mcStat:
                        err *= pulls["prop_binBin0_bin"+str(i_bin)].sigma
#                        err = err**pulls["prop_binBin0_bin"+str(i_bin)].sigma
                    shapeH.SetBinContent( i_bin+1, err )
                    uncertainties[_bin][est][unc] = err
                if not "histo" in uncertainties.keys(): uncertainties["histo"] = {}
                if not est in uncertainties["histo"].keys(): uncertainties["histo"][est] = {}
                uncertainties["histo"][est][unc] = shapeH.Clone()
            else:
                shapeH.Add(shapes[est],-1)
                shapeH.Divide(shapes[est])
                if postFit:
                    shapeH.Scale(pulls[unc].sigma)
#                    for i_bin in range(shapeH.GetNbinsX()):
#                        shapeH.SetBinContent(i_bin+1, shapeH.GetBinContent(i_bin+1)**pulls[unc].sigma)
                for i_bin in range(shapeH.GetNbinsX()):
                    _bin = "Bin"+str(i_bin)
                    if not _bin in uncertainties.keys(): uncertainties[_bin] = {}
                    if not est in uncertainties[_bin].keys(): uncertainties[_bin][est] = {}
                    uncertainties[_bin][est][unc] = shapeH.GetBinContent( i_bin+1 )
                if not "histo" in uncertainties.keys(): uncertainties["histo"] = {}
                if not est in uncertainties["histo"].keys(): uncertainties["histo"][est] = {}
                uncertainties["histo"][est][unc] = shapeH.Clone()
        
        self.uncertainties[key] = uncertainties

        if bin or estimate or nuisance or systOnly:
            return self.__filterDict( self.uncertainties[key], bin=bin, estimate=estimate, nuisance=nuisance, systOnly=systOnly )
        else:
            return self.uncertainties[key]


    def getUncertaintiesFromTxtCard( self, bin=None, estimate=None, nuisance=None, postFit=False, systOnly=False ):
        # uncertainties from txt card
        # return safed uncertainties if available
        key = "postFit" if postFit else "preFit"
        if self.uncertainties[key]:
            if bin or estimate or nuisance or systOnly:
                return self.__filterDict( self.uncertainties[key], bin=bin, estimate=estimate, nuisance=nuisance, systOnly=systOnly )
            else:
                return self.uncertainties[key]

        allUnc        = self.getNuisancesList( systOnly=systOnly )
        allEst        = self.getProcessList( unique=False )
        rateParams    = self.getRateParameter()
        pulls         = self.getPulls( postFit=postFit )
        binList       = self.getBinList( unique=False )
        uncertainties = {}

        with open( self.txtCard ) as f:
            for line in f:
                if not line.split() or line.startswith("-"): continue
                unc = line.split()[0] 
                if unc not in allUnc: continue
                if unc in rateParams.keys(): continue # remove rate parameters as they would be 0 anyway
                for i_bin, _bin in enumerate(binList):
                    est = allEst[i_bin]
                    if not _bin in uncertainties.keys(): uncertainties[_bin] = {}
                    if not est in uncertainties[_bin].keys(): uncertainties[_bin][est] = {}
                    try:
                        uncertainties[_bin][est][unc] = float(line.split()[2:][i_bin])-1
                    except:
                        uncertainties[_bin][est][unc] = 0
                    if postFit and uncertainties[_bin][est][unc]:
                        uncertainties[_bin][est][unc] *= pulls[unc].sigma
#                        uncertainties[_bin][est][unc] = uncertainties[_bin][est][unc]**pulls[unc].sigma

        self.uncertainties[key] = uncertainties

        if bin or estimate or nuisance or systOnly:
            return self.__filterDict( self.uncertainties[key], bin=bin, estimate=estimate, nuisance=nuisance, systOnly=systOnly )
        else:
            return self.uncertainties[key]

    def getObservation( self, bin=None, directory="total" ):
        return {dir:{ b:b_dict["data"] for b, b_dict in o.iteritems() } for dir, o in self.getEstimates( postFit=False, bin=bin, estimate="data", directory=directory ).iteritems()}

    def getEstimates( self, bin=None, estimate=None, directory="total", postFit=False ):
        key    = "postFit" if postFit else "preFit"
        if self.estimates[key]:
            ests = self.estimates[key]
            all = { d:self.__filterDict( dic, bin=bin, estimate=estimate ) if bin else dic for d, dic in ests.iteritems() } 
            if directory: return {directory:all[directory]}                
            else: return all

        regionHistos = self.getRegionHistos( postFit=postFit, plotBins=None, directory=None, addRateUncertainty=False )
        processes    = self.getProcessesPerBin( bin=None )
        yields       = {}
        for dir, histoDict in regionHistos.iteritems():
            tmp          = {}
            yields[dir]  = {}
            for est, h in histoDict.iteritems():
#                if est not in processes and est != "data" and not "total" in est: continue
                tmp[est] = {}
                for i in range(h.GetNbinsX()):
                    y = h.GetBinContent(i+1)
                    e = h.GetBinError(i+1) # Attention, this is the full error, not only MC statistics! Fixme
                    key = dir + "_Bin%i"%i if directory != "total" else "Bin0_Bin%i"%i
                    if key in tmp[est].keys(): tmp[est][key] += u_float( y, e )
                    else:                      tmp[est][key]  = u_float( y, e )
                    yields[dir][key]    = {}

            # stupid restructuring to make it compatible w/ other functions
            for b in yields[dir].keys():
                for est in tmp.keys():
                    yields[dir][b][est] = tmp[est][b]

        self.estimates[key] = yields

        all = { d:self.__filterDict( dic, bin=bin, estimate=estimate ) if bin else dic for d, dic in yields.iteritems() } 
        if directory: return {directory:all[directory]}                
        else: return all

    def getNuisanceHistosFromShapeCard( self, postFit=False, plotBins=None, nuisances=None, directory=None ):

        if any( ["prop" in n for n in nuisances] ):
            nuisances = [ n for n in nuisances if not "prop" in n ] + ["stat"]

        allEst   = self.getProcessList( unique=True )
        histDict = self.getUncertaintiesFromShapeCard( postFit=postFit )["histo"]
        yields   = self.__getShapeObject( key=None )
        rateParams = self.getRateParameter( postFit=True ).keys()
        nuisanceHistos = {}

        dirs = self.channels if not directory else [directory]

        for dir in dirs:
            regions = self.getRegionHistos( postFit=postFit, plotBins=plotBins, directory=dir, addRateUncertainty=False )[dir]
            nuisanceHistos[dir] = {}
            for i_n, nuisance in enumerate(nuisances):

                if nuisance in rateParams or "Lumi" in nuisance: #quick fix, use getNuisanceHistos for lnN, FIXME
                    nuisanceHistos[dir][nuisance] = self.getNuisanceHistos( postFit=postFit, plotBins=plotBins, nuisances=[nuisance], directory=directory )[dir][nuisance]
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
                nuisanceHistos[dir][nuisance]["up"].style        = styles.lineStyle( ROOT.kSpring-1-i_n, width=3 ) #change to dynamic style
                nuisanceHistos[dir][nuisance]["down"].style      = styles.lineStyle( ROOT.kOrange+7+i_n, width=3 )
                nuisanceHistos[dir][nuisance]["up"].legendText   = nuisance + " (+1#sigma)"
                nuisanceHistos[dir][nuisance]["down"].legendText = nuisance + " (-1#sigma)"

        return nuisanceHistos

    def getNuisanceHistos( self, postFit=False, plotBins=None, nuisances=None, directory="total" ):

        histDict       = self.getRegionHistos( postFit=postFit, plotBins=None, directory=directory, addRateUncertainty=False )
        nuisanceHistos = {}

        if directory == "total":
            binList        = self.getBinList( unique=True )

        for dir, hist in histDict.iteritems():
            nuisanceHistUp      = hist["total_signal"].Clone()
            nuisanceHistDown    = hist["total_signal"].Clone()
            nuisanceHistos[dir] = {}
            for i_n, nuisance in enumerate(nuisances):
                nuisanceYields  = self.getNuisanceYields( nuisance, postFit=postFit, directory=directory )
                for i in range(nuisanceHistUp.GetNbinsX()):
                    if self.year == "combined" and directory != "total":
                        key = directory+"_Bin"+str(i)
                    elif self.year == "combined" and directory == "total":
                        key = binList[i]
                    else:
                        key = "Bin"+str(i)

                    nDict = nuisanceYields[key]
                    nuisanceHistUp.SetBinContent(   i+1, nDict["up"] )
                    nuisanceHistDown.SetBinContent( i+1, nDict["down"] )
    
                nuisanceHistos[dir][nuisance]                    = { "up":nuisanceHistUp.Clone(), "down":nuisanceHistDown.Clone() }
                nuisanceHistos[dir][nuisance]["up"].style        = styles.lineStyle( ROOT.kSpring-1-i_n, width=3 ) #change to dynamic style
                nuisanceHistos[dir][nuisance]["down"].style      = styles.lineStyle( ROOT.kOrange+7+i_n, width=3 )
                nuisanceHistos[dir][nuisance]["up"].legendText   = nuisance + " (+1#sigma)"
                nuisanceHistos[dir][nuisance]["down"].legendText = nuisance + " (-1#sigma)"

        return nuisanceHistos

    def createRebinnedResults( self, rebinningCardFile, postfit=False ):

        # create environment
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp/", ustr)
        print "Creating "+uniqueDirname
        os.makedirs(uniqueDirname)

        resPath      = self.txtCard.replace(".txt","/")
        if not os.path.isdir(resPath): os.makedirs(resPath)
        if self.year == "combined":
            if not os.path.isdir(resPath.replace("COMBINED","2016")): os.makedirs(resPath.replace("COMBINED","2016"))
            if not os.path.isdir(resPath.replace("COMBINED","2017")): os.makedirs(resPath.replace("COMBINED","2017"))
            if not os.path.isdir(resPath.replace("COMBINED","2018")): os.makedirs(resPath.replace("COMBINED","2018"))

        f            = "masked_" + rebinningCardFile.split("/")[-1]
        resTxtFile   = os.path.join( resPath, f )
        resShapeFile = os.path.join( resPath, f.replace(".txt","_shapeCard.txt") )
        resShapeRoot = os.path.join( resPath, f.replace(".txt","_shape.root") )
        shutil.copyfile(os.path.join(os.environ['CMSSW_BASE'], 'src', 'Analysis', 'Tools', 'python', 'cardFileWriter', 'diffNuisances.py'), os.path.join(uniqueDirname, 'diffNuisances.py'))

        # combine fit card and muted card
        print "combining cards for muted fit"
        if self.year == "combined":
            cmd  = "cd "+uniqueDirname+";combineCards.py fit_dc_2016=%s fit_dc_2017=%s fit_dc_2018=%s dc_2016=%s dc_2017=%s dc_2018=%s > combinedCard.txt; text2workspace.py combinedCard.txt --channel-masks"%(self.shapeCard.replace("COMBINED","2016"), self.shapeCard.replace("COMBINED","2017"), self.shapeCard.replace("COMBINED","2018"), rebinningCardFile.replace(".txt","_shapeCard.txt").replace("COMBINED","2016"), rebinningCardFile.replace(".txt","_shapeCard.txt").replace("COMBINED","2017"), rebinningCardFile.replace(".txt","_shapeCard.txt").replace("COMBINED","2018"))
        else:
            cmd  = "cd "+uniqueDirname+";combineCards.py fit=%s Bin0=%s > combinedCard.txt; text2workspace.py combinedCard.txt --channel-masks"%(self.shapeCard, rebinningCardFile.replace(".txt","_shapeCard.txt"))
        print "Executing command: %s"%cmd
        os.system(cmd)

        # combine fit card and muted card, create workspace
        print "run text2workspace"
        cmd  = "cd "+uniqueDirname+";text2workspace.py combinedCard.txt --channel-masks --X-allow-no-signal -m 125"
        print "Executing command: %s"%cmd
        os.system(cmd)

        # also combine txt cards
        print "combining txt cards for muted fit"
        if self.year == "combined":
            cmd  = "cd "+uniqueDirname+";combineCards.py fit_dc_2016=%s fit_dc_2017=%s fit_dc_2018=%s dc_2016=%s dc_2017=%s dc_2018=%s > txtCard.txt"%(self.txtCard.replace("COMBINED","2016"), self.txtCard.replace("COMBINED","2017"), self.txtCard.replace("COMBINED","2018"), rebinningCardFile.replace("COMBINED","2016"), rebinningCardFile.replace("COMBINED","2017"), rebinningCardFile.replace("COMBINED","2018"))
        else:
            cmd  = "cd "+uniqueDirname+";combineCards.py fit=%s Bin0=%s > txtCard.txt"%(self.txtCard, rebinningCardFile)
        print "Executing command: %s"%cmd
        os.system(cmd)

        # run fit with masked (muted) card
        print "run FitDiagnostics"
        if self.year == "combined":
            cmd = "cd "+uniqueDirname+";combine combinedCard.root --robustHesse 1 --forceRecreateNLL -M FitDiagnostics --saveShapes --saveNormalizations --saveOverall --saveWithUncertainties --setParameters mask_dc_2016=1,mask_dc_2017=1,mask_dc_2018=1"
        else:
            cmd = "cd "+uniqueDirname+";combine combinedCard.root --robustHesse 1 --forceRecreateNLL -M FitDiagnostics --saveShapes --saveNormalizations --saveOverall --saveWithUncertainties --setParameters mask_Bin0=1"
        print "Executing command: %s"%cmd
        os.system(cmd)

        # copy cards to final location
        logger.info("Putting combined card into dir %s", resPath)
        if self.year == "combined":
            shutil.copyfile(rebinningCardFile.replace("COMBINED","2016"),                                  resTxtFile.replace("COMBINED","2016"))
            shutil.copyfile(rebinningCardFile.replace(".txt","_shapeCard.txt").replace("COMBINED","2016"), resShapeFile.replace("COMBINED","2016"))
            shutil.copyfile(rebinningCardFile.replace(".txt","_shape.root").replace("COMBINED","2016"),    resShapeRoot.replace("COMBINED","2016"))
            shutil.copyfile(rebinningCardFile.replace("COMBINED","2017"),                                  resTxtFile.replace("COMBINED","2017"))
            shutil.copyfile(rebinningCardFile.replace(".txt","_shapeCard.txt").replace("COMBINED","2017"), resShapeFile.replace("COMBINED","2017"))
            shutil.copyfile(rebinningCardFile.replace(".txt","_shape.root").replace("COMBINED","2017"),    resShapeRoot.replace("COMBINED","2017"))
            shutil.copyfile(rebinningCardFile.replace("COMBINED","2018"),                                  resTxtFile.replace("COMBINED","2018"))
            shutil.copyfile(rebinningCardFile.replace(".txt","_shapeCard.txt").replace("COMBINED","2018"), resShapeFile.replace("COMBINED","2018"))
            shutil.copyfile(rebinningCardFile.replace(".txt","_shape.root").replace("COMBINED","2018"),    resShapeRoot.replace("COMBINED","2018"))

        shutil.copyfile(rebinningCardFile,                                  resTxtFile)
        shutil.copyfile(rebinningCardFile.replace(".txt","_shapeCard.txt"), resShapeFile)
        shutil.copyfile(rebinningCardFile.replace(".txt","_shape.root"),    resShapeRoot)
        shutil.copyfile(uniqueDirname+"/txtCard.txt",         resTxtFile.replace(  "masked","masked_fit"))
        shutil.copyfile(uniqueDirname+"/combinedCard.txt",    resShapeFile.replace("masked","masked_fit"))
        shutil.copyfile(uniqueDirname+"/combinedCard.root",   resShapeFile.replace(".txt",".root"))
        shutil.copyfile(uniqueDirname+"/fitDiagnostics.root", resShapeFile.replace(".txt","_FD.root"))

        os.remove(uniqueDirname+"/fitDiagnostics.root")

        # get number of bins of the rebinning card
        rbResults = CombineResults( cardFile=rebinningCardFile, plotDirectory=self.plotDirectory, year=self.year, bkgOnly=self.bkgOnly, isSearch=self.isSearch )
        newPulls  = rbResults.getPulls().keys()
        nBins = len( rbResults.getBinList( unique=True ) )
        del rbResults

        # run fit with masked (muted) card
        print "run FitDiagnostics stat only"
        if self.year == "combined":
            cmd = "cd "+uniqueDirname+";combine combinedCard.root --profilingMode none -M FitDiagnostics --saveWithUncertainties --saveShapes --saveNormalizations --saveOverall --setParameters mask_dc_2016=1,mask_dc_2017=1,mask_dc_2018=1"
        else:
            cmd = "cd "+uniqueDirname+";combine combinedCard.root --profilingMode none -M FitDiagnostics --saveWithUncertainties --saveShapes --saveNormalizations --saveOverall --setParameters mask_Bin0=1"
        print "Executing command: %s"%cmd
        os.system(cmd)

        shutil.copyfile(uniqueDirname+"/fitDiagnostics.root", resShapeFile.replace(".txt","_statOnly_FD.root"))

        shutil.rmtree(uniqueDirname)

        # rewrite content in a similar way to the combine fit results
        self.__rewriteRebinnedFile( resShapeFile.replace(".txt","_FD.root"), postfit=postfit, nBins=nBins )
        self.__rewriteRebinnedFile( resShapeFile.replace(".txt","_statOnly_FD.root"), postfit=postfit, nBins=nBins, statOnly=True )

        return resTxtFile

    def getRegionHistos( self, postFit=False, plotBins=None, nuisances=None, bkgSubstracted=False, labelFormater=None, directory=None, addStatOnlyHistos=False, addRateUncertainty=True ):
        hists = self.__regionHistos( postFit=postFit, plotBins=plotBins, nuisances=nuisances, bkgSubstracted=bkgSubstracted, labelFormater=labelFormater, statOnly=False, addRateUncertainty=addRateUncertainty )
        if addStatOnlyHistos:
            hists_stat = self.__regionHistos( postFit=postFit, plotBins=plotBins, nuisances=nuisances, bkgSubstracted=bkgSubstracted, labelFormater=labelFormater, statOnly=True,  addRateUncertainty=False )
            for dir, dic in hists_stat.iteritems():
                for h_key, h in dic.iteritems():
                    hists[dir][h_key+"_stat"] = h.Clone()

        return {directory:hists[directory]} if directory else hists

    def __regionHistos( self, postFit=False, plotBins=None, nuisances=None, bkgSubstracted=False, labelFormater=None, statOnly=False, addRateUncertainty=True ):

        hists    = {}
        key    = "postFit" if postFit else "preFit"
        subkey = self.__getSubKey( plotBins=plotBins )

        if not statOnly and subkey in self.regionHistos[key].keys() and self.regionHistos[key][subkey]:
            hists = self.regionHistos[key][subkey]
            return hists

        if   postFit and not self.bkgOnly: dirName = "shapes_fit_s"
        elif postFit and     self.bkgOnly: dirName = "shapes_fit_b"
        else:                              dirName = "shapes_prefit"

        if statOnly: fit = self.__getStatOnlyFitObject( key=dirName )
        else:        fit = self.__getFitObject( key=dirName )
        histDict = { ch: [ x.GetName() for x in fit.Get(ch).GetListOfKeys() if x.GetName() != "data" ] + [ "data" ] for ch in self.channels }

        for dir, histList in histDict.iteritems():
            hists[dir]    = {}
            histList.sort()
            histList = filter( lambda hist: "total_covar" not in hist and "process_" not in hist, histList )
            for hist in histList:

                hists[dir][hist] = fit.Get(dir+"/"+hist).Clone()

                # change TGraph type to TH1F type for data
                if "data" in hist:# and not "_rebinned" in self.txtCard:
                    dataHist = hists[dir][histList[0]].Clone()
                    dataHist.Reset()
                    dataHist.SetName("data")

                    if type( hists[dir][hist] ) == ROOT.TGraphAsymmErrors:
                        for i in range(dataHist.GetNbinsX()):
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

                if self.year == "combined":
                    k = "data" if "data" in hist else hist
                    hists[dir][k].GetXaxis().SetRangeUser(0, int(hists[dir][k].GetNbinsX()/3.))

            if nuisances and not bkgSubstracted: # currently no single-nuisance plots with bkg substracted histograms #FIXME
                if isinstance( nuisances, str ): nuisances = [nuisances]
                hists[dir].update( self.getNuisanceHistosFromShapeCard( postFit=postFit, plotBins=None, nuisances=nuisances, directory=dir )[dir] )

            # add rate parameter postfit uncertainty to total histogram and total_background histogram as it is not included in combine
            if postFit and addRateUncertainty and not statOnly:
                tot        = "total" if "total" in hists[dir].keys() else "total_overall"
                rateParams = self.getRateParameter( postFit=True ).keys()
                totalHist  = hists[dir][tot]
                for rateParam in rateParams:
                    rateHisto = self.getNuisanceHistosFromShapeCard( postFit=True, nuisances=[rateParam], directory=dir )[dir][rateParam]["up"].Clone()
                    rateHisto.Add( totalHist, -1 )
                    for i in range(hists[dir][tot].GetNbinsX()):
                        err  = rateHisto.GetBinContent(i+1)
                        terr = hists[dir][tot].GetBinError(i+1)
                        hists[dir][tot].SetBinError(i+1, math.sqrt(err*err+terr*terr))
                        terr = hists[dir]["total_background"].GetBinError(i+1)
                        hists[dir]["total_background"].SetBinError(i+1, math.sqrt(err*err+terr*terr))
       
            labels = self.getBinLabels( labelFormater=labelFormater )
            if labels:
                for h_key, h in hists[dir].iteritems():
                    if isinstance( h, dict ):
                        for i in range(h["up"].GetNbinsX()):
                            h["up"].GetXaxis().SetBinLabel( i+1, labels[i] )
                            h["down"].GetXaxis().SetBinLabel( i+1, labels[i] )
                        h["up"].LabelsOption("v","X") #"vu" for 45 degree labels
                        h["down"].LabelsOption("v","X") #"vu" for 45 degree labels
                    else:
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

        if not statOnly: self.regionHistos[key][subkey] = hists

        if bkgSubstracted:
            for dir, histList in histDict.iteritems():
                tot = "total" if "total" in hists[dir].keys() else "total_overall"

                # remove error on total background
                for b in range(hists[dir]["total_background"].GetNbinsX()):
                    hists[dir]["total_background"].SetBinError(b+1, 0)

                # use total - bkg as signal to get the full uncertainty
                hists[dir] = {"data":hists[dir]["data"], "signal":hists[dir][tot].Clone(),"total":hists[dir][tot].Clone(),"total_background":hists[dir]["total_background"]}

                for i in range( hists[dir]["signal"].GetNbinsX() ):
                    print statOnly, i, hists[dir]["signal"].GetBinError(i+1) / hists[dir]["signal"].GetBinContent(i+1), hists[dir][tot].GetBinError(i+1) / hists[dir][tot].GetBinContent(i+1)
                print
                print
                print
                hists[dir]["data"].Add( hists[dir]["total_background"], -1 )
                hists[dir]["signal"].Add( hists[dir]["total_background"], -1 )
                hists[dir]["total_background"].Scale(0)
                for i in range( hists[dir]["signal"].GetNbinsX() ):
                    print statOnly, i, hists[dir]["signal"].GetBinError(i+1) / hists[dir]["signal"].GetBinContent(i+1), hists[dir][tot].GetBinError(i+1) / hists[dir][tot].GetBinContent(i+1)
        return hists

    def getRegionHistoList( self, regionHistos, processes=None, noData=False, sorted=False, bkgSubstracted=False, directory="total"):
        # get the list of histograms and the ratio list for plotting a region plot

        if bkgSubstracted: return [ [regionHistos["signal"]], [regionHistos["data"]] ], [(1,0),(0,0)]

        for p in processes:
            if not p in regionHistos.keys():
                # some histograms are 0, still should be in the legend
                logger.info("Histogram for %s not found! Creating one and setting it to 0! Continuing..."%p)
                regionHistos[p] = regionHistos["signal"].Clone()
                self.__copyHistoSettings( fromHist=regionHistos["signal"], toHist=regionHistos[p], plotBins=None )
                regionHistos[p].Scale(0.)
                regionHistos[p].SetName(p)
                del regionHistos[p].legendText
                regionHistos[p].notInLegend = True
    

        nuisances    = self.getNuisancesList( systOnly=False ) + ["totalSyst"]
        binProcesses = self.getProcessesPerBin()
        ratioHistos  = []
        bins         = len(self.getBinLabels()) if self.year != "combined" else int(len(self.getBinLabels())/3.)
        i_n          = 0

        if sorted:
            histoList = [[]]
            for i in range( bins ):
                proc_list = []
                key = directory + "_Bin%i"%i 
                if directory == "total" or "Bin0":
                    key = "Bin%i"%i
                elif directory != "total" and any( ["Bin" in k for k in binProcesses.keys()]):
                    key = directory + "_Bin%i"%i
                else:
                    key = directory
                for p in binProcesses[key]:

                    if p in regionHistos.keys():
                        # set only one bin != 0
                        tmp = regionHistos[p].Clone( p + "_Bin%i_%s"%(i, str(uuid.uuid4())) )
                        tmp.Scale(0.)
                        tmp.SetBinContent( i+1, regionHistos[p].GetBinContent(i+1) )
                        self.__copyHistoSettings( fromHist=regionHistos[p], toHist=tmp, plotBins=None )
                    else:
                        tmp = regionHistos["signal"].Clone( p + "_Bin%i_%s"%(i, str(uuid.uuid4())) )
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
#                histoList[0] += copy.copy(proc_list)
                histoList[0] += proc_list

        else:
            histoList = [ [p_h for p, p_h in regionHistos.iteritems() if p in binProcesses["Bin0"] ] ]
            histoList[0].sort( key=lambda h: -regionHistos[p].Integral() )

        # add data histos
        if not noData:
            histoList   += [ [regionHistos["data"]] ]
            ratioHistos += [ (1,0) ]
            i_n         += 1

        # add nuisance histos at last
        for n in nuisances:
            if n in regionHistos.keys() and isinstance( regionHistos[n], dict ):
                histoList   += [ [regionHistos[n]["up"]], [regionHistos[n]["down"]] ]
                ratioHistos += [ ((i_n)*2,0),((i_n)*2+1,0) ]
                i_n         += 1

        for i in range( regionHistos["signal"].GetNbinsX() ):
            for h_list in histoList:
                for h in h_list:
                    # make that more dynamic FIXME
                    h.GetXaxis().SetBinLabel( i+1, regionHistos["signal"].GetXaxis().GetBinLabel( i+1 ) )
                    h.LabelsOption("v","X") #"vu" for 45 degree labels

        return histoList, ratioHistos
        

    def setPlotDirectory( self, plotDirectory ):
        self.plotDirectory = plotDirectory

    def getImpactPlot( self, expected=False, printPNG=False, cores=1 ):

        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join("/tmp", ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

        plotName = "impacts"
        if self.bkgOnly: plotName += "_bkgOnly"
        if expected:     plotName += "_expected"

        if self.bkgOnly: options = "--freezeParameters r --setParameters r=%i"%(0 if self.isSearch else 1)
        else:            options = "--rMin 0 --rMax 10"

        cd             = "cd %s"%uniqueDirname
        robustFit      = "combineTool.py -M Impacts -m 125 -d %s --doInitialFit --robustFit 1 %s"%(self.rootWorkSpace, options)
        impactFits     = "combineTool.py -M Impacts -m 125 -d %s --robustFit 1 --doFits --parallel %i %s"%( self.rootWorkSpace, cores, options )
        extractImpact  = "combineTool.py -M Impacts -m 125 -d %s -o impacts.json"%self.rootWorkSpace
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

    def getCorrelationHisto( self, systOnly=False ):

        if not self.fitResult:
            raise ValueError( "Root file of fit result not found! Running in limited mode, thus cannot get the object needed!" )

        if self.correlationHisto:
            if systOnly:
                nuisSyst = self.getNuisancesList( systOnly=True )
                nuisAll  = self.getNuisancesList( systOnly=False )
                hist     = self.correlationHisto.Clone(str(uuid.uuid4()))
                hist.GetXaxis().SetRangeUser(0,len(nuisSyst))
                hist.GetYaxis().SetRangeUser(len(nuisAll)-len(nuisSyst),len(nuisAll))
                hist.LabelsOption("v","X")
                return hist
            return self.correlationHisto

        fit      = self.__getFitObject( key="fit_b" if self.bkgOnly else "fit_s" )
        corrhist = copy.deepcopy(fit.correlationHist())

        # bit of formating
        corrhist.GetZaxis().SetRangeUser(-1,1)
        corrhist.LabelsOption("v","X")

        self.correlationHisto = corrhist

        if systOnly:
            nuisSyst = self.getNuisancesList( systOnly=True )
            nuisAll  = self.getNuisancesList( systOnly=False )
            hist     = corrhist.Clone(str(uuid.uuid4()))
            hist.GetXaxis().SetRangeUser(0,len(nuisSyst))
            hist.GetYaxis().SetRangeUser(len(nuisAll)-len(nuisSyst)+1,len(nuisAll)+1)
            return hist
        return self.correlationHisto

    def getCovarianceHisto( self, labelFormater=None, postFit=False ):
        # get the TH2D covariance matrix plot

        if not self.fitResult:
            raise ValueError( "Root file of fit result not found! Running in limited mode, thus cannot get the object needed!" )

        key = "postFit" if postFit else "preFit"
        if self.covarianceHistos[key]:
            return self.covarianceHistos[key]

        if postFit:
            dirName = "shapes_fit_b" if self.bkgOnly else "shapes_fit_s"
        else:
            dirName = "shapes_prefit"

        fit = self.__getFitObject( key=dirName )
        self.covarianceHistos[key] = copy.deepcopy( fit.Get("overall_total_covar") )
        # normalize
        self.covarianceHistos[key].Scale(1./self.covarianceHistos[key].GetMaximum())

        # set labels
        self.covarianceHistos[key].LabelsOption("v","X")
        labels = self.getBinLabels( labelFormater=labelFormater )
        for i in range(self.covarianceHistos[key].GetNbinsY()):
            self.covarianceHistos[key].GetYaxis().SetBinLabel( i+1, labels[i] )
            self.covarianceHistos[key].GetXaxis().SetBinLabel( i+1, labels[i] )

        return self.covarianceHistos[key]

