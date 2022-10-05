import ROOT

from array import *

import shutil
import os
from math import log

from Analysis.Tools.helpers import writeObjToFile

# Logging
import logging
logger = logging.getLogger(__name__)

import re
def natural_sort(list, key=lambda s:s):
    """
    Sort the list into natural alphanumeric order.
    http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
    """
    def get_alphanum_key_func(key):
        convert = lambda text: int(text) if text.isdigit() else text
        return lambda s: [convert(c) for c in re.split('([0-9]+)', key(s))]
    sort_key = get_alphanum_key_func(key)

    lc = sorted(list, key=sort_key)
    return lc

class cardFileWriter:
    def __init__(self):
        self.reset()
        self.releaseLocation = os.path.abspath('.')

    def reset(self):
        self.bins = []
        self.muted = {}
        self.uncertainties = []
        self.uncertaintyVal = {}
        self.uncertaintyString = {}
        self.processes = {}
        self.expectation = {}
        self.observation = {}
        self.contamination = {}
        self.niceNames = {}
        self.defWidth = 12
        self.precision = 7
        self.maxUncNameWidth = 30
        self.maxUncStrWidth = 30
        self.hasContamination=False
        self.rateParameters = []
        self.freeParameters = []
        self.precision = 10
	self.freeParamBins={}
	self.regionMapping = {}
	self.CR = []
	self.SR = []
    def addFreeParamBin(self,name,binlist):
	self.freeParamBins= {name:binlist}

    def addMap(self,regionMapping):
	self.regionMapping=regionMapping
    def addCR(self,CRbins):
        self.CR=CRbins
        self.nCR = len(CRbins)

    def addSR(self,SRbins):
	self.SR=SRbins

    def setPrecision(self, prec):
        self.precision = prec

    def addBin(self, name, processes, niceName=""):
        if len(name)>30:
            logger.info("Name for bin",name,"too long. Max. length is 30.")
            return
        if self.niceNames.has_key(name):
            logger.info("Bin already there! (",name,")")
            return
        for p in processes:
            if len(p)>30:
                logger.info("Name for process", p, "in bin", name, "is too long. Max. length is 30.")
                return
        self.niceNames[name]=niceName
        self.bins.append(name)
        self.muted[name]=False
        self.processes[name] = ["signal"]+processes

    def addUncertainty(self, name, t, n=0):
        assert len(name)<self.maxUncNameWidth,  "That's too long: %s. Max. length is %i"%(name, self.maxUncNameWidth)
        if self.uncertainties.count(name):
            print "Uncertainty already there! (",name,")"
            return
        self.uncertainties.append(name)
        self.uncertaintyString[name] = t
        if t=="gmN":
            if n==0:
                print "gmN Uncertainty with n=0! Specify n as third argument: addUncertainty(..., 'gmN', n)"
                return
            self.uncertaintyString[name] = t+' '+str(n)
        if len(self.uncertaintyString[name])>self.maxUncStrWidth:
            print "That's too long:",self.uncertaintyString[name],"Max. length is", self.maxUncStrWidth
            del self.uncertaintyString[name]
            return
    
    def addRateParameter(self, p, value, r):
        if [ a[0] for a in self.rateParameters ].count(p):
            logger.info("Rate parameter for process %s already added!"%p)
            return
        self.rateParameters.append((p, value, r))

    def addFreeParameter(self, name, p, value, r):
        if [ a[0] for a in self.freeParameters ].count(name):
            logger.info("Free parameter for process %s already added!"%p)
            return
        self.freeParameters.append((name, p, value, r))

    def specifyExpectation(self, b, p, exp):
        self.expectation[(b,p)] = round(exp, self.precision)

    def specifyObservation(self, b, obs):
        if not isinstance(obs, int):
            print "Observation not an integer! (",obs,")"
            return
        self.observation[b] = obs

    def specifyContamination(self, b, cont):
        self.contamination[b] = cont
        self.hasContamination = True

    def specifyFlatUncertainty(self, u,  val):
        if u not in self.uncertainties:
            print "This uncertainty has not been added yet!",u,"Available:",self.uncertainties
            return
        print "Adding ",u,"=",val,"for all bins and processes!"
        for b in self.bins:
            for p in self.processes[b]:
                self.uncertaintyVal[(u,b,p)] = round(val,self.precision)

    def specifyUncertainty(self, u, b, p, val):
        if u not in self.uncertainties:
            print "This uncertainty has not been added yet!",u,"Available:",self.uncertainties
            return
        if b not in self.bins:
            print "This bin has not been added yet!",b,"Available:",self.bins
            return
        if p not in self.processes[b]:
            print "Process ", p," is not in bin",b,". Available for ", b,":",self.processes[b]
            return
        if val<0:
            print "Warning! Found negative uncertainty %f for yield %f in %r. Reversing sign under the assumption that the correlation pattern is irrelevant (check!)."%(val, self.expectation[(b, p)], (u,b,p))
            _val=1.0
        else:
            _val = val
        self.uncertaintyVal[(u,b,p)] = round(_val,self.precision)

    def getUncertaintyString(self, k):
        u, b, p = k
        if self.uncertaintyString[u].count('gmN'):
            if self.uncertaintyVal.has_key((u,b,p)) and self.uncertaintyVal[(u,b,p)]>0.:
                n = float(self.uncertaintyString[u].split(" ")[1])
                return self.mfs(self.expectation[(b, p)]/float(n))
            else: return '-'
        if self.uncertaintyVal.has_key((u,b,p)):
            return self.mfs(self.uncertaintyVal[(u,b,p)])
        return '-'

    def checkCompleteness(self):
        for b in self.bins:
            if not self.observation.has_key(b) or not self.observation[b]<float('inf'):
                print "No valid observation for bin",b
                return False
            if self.hasContamination and (not self.contamination.has_key(b) or not self.contamination[b] < float('inf')):
                print "No valid contamination for bin",b
                return False
            if len(self.processes[b])==0:
                print "Warning, bin",b,"has no processes"
            for p in self.processes[b]:
                if not self.expectation.has_key((b,p)) or not self.expectation[(b,p)]<float('inf'):
                    print "No valid expectation for bin/process ",(b,p)
                    return False
            for k in self.uncertaintyVal.keys():
                if not self.uncertaintyVal[k]<float('inf'):
                    print "Uncertainty invalid for",k,':',self.uncertaintyVal[k]
                    return False
        return True

    def mfs(self, f):
        return str(round(float(f),self.precision))

    def writeToFile(self, fname, shapeFile=False, noMCStat=False):

        print "noMCStat: {}".format(noMCStat)

        import datetime, os
        if not self.checkCompleteness():
            print "Incomplete specification."
            return
        allProcesses=[]
        numberID = {}
        i=1
        for b in self.bins:
            if not self.muted[b]:
                for p in self.processes[b]:
                    if not p in allProcesses and not p=='signal':
                        allProcesses.append(p)
                        numberID[p] = i
                        i+=1
        unmutedBins = [b for b in self.bins if not self.muted[b]]
        nBins = len(unmutedBins)
        numberID['signal'] = 0
        lspace = (self.maxUncStrWidth + self.maxUncNameWidth + 2)
        if not os.path.exists(os.path.dirname(fname)):
            os.makedirs(os.path.dirname(fname))
        outfile = file(fname, 'w')
        outfile.write('#cardFileWriter, %s'%datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")+'\n')
        outfile.write('imax %i'%nBins+'\n')
        outfile.write('jmax *\n')
        outfile.write('kmax *\n')
        outfile.write('\n')

        for b in self.bins:
            if not self.muted[b]:
                outfile.write( '# '+b+': '+self.niceNames[b]+'\n')
            else:
                outfile.write( '#Muted: '+b+': '+self.niceNames[b]+'\n')
        outfile.write( '\n')

        if shapeFile:
            outfile.write( 'shapes * * %s  $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC \n'%shapeFile)
            outfile.write( '\n')

        outfile.write( 'bin'.ljust(lspace)              +(' '.join([b.rjust(self.defWidth) for b in unmutedBins] ) ) +'\n')
        outfile.write( 'observation'.ljust(lspace)      +(' '.join([str(self.observation[b]).rjust(self.defWidth) for b in unmutedBins]) )+'\n')
        if self.hasContamination:
            outfile.write( 'contamination'.ljust(lspace)  +(' '.join([str(self.contamination[b]).rjust(self.defWidth) for b in unmutedBins]) )+'\n')
        outfile.write('\n')
        outfile.write( 'bin'.ljust(lspace)              +(' '.join( [' '.join([b.rjust(self.defWidth) for p in self.processes[b]] ) for b in unmutedBins]) ) +'\n')
        outfile.write( 'process'.ljust(lspace)          +(' '.join( [' '.join([p.rjust(self.defWidth) for p in self.processes[b]] ) for b in unmutedBins]) ) +'\n')
        outfile.write( 'process'.ljust(lspace)          +(' '.join( [' '.join([str(numberID[p]).rjust(self.defWidth) for p in self.processes[b]] ) for b in unmutedBins]) ) +'\n')
        outfile.write( 'rate'.ljust(lspace)             +(' '.join( [' '.join([self.mfs(self.expectation[(b,p)]).rjust(self.defWidth) for p in self.processes[b]] ) for b in unmutedBins]) )+'\n')
        outfile.write('\n')

        for u in self.uncertainties:
            outfile.write( u.ljust(self.maxUncNameWidth)+' '+self.uncertaintyString[u].ljust(self.maxUncStrWidth)+' '+
                                          ' '.join( [' '.join([self.getUncertaintyString((u,b,p)).rjust(self.defWidth) for p in self.processes[b]] ) for b in unmutedBins]) +'\n')

	print self.rateParameters
        for p in self.rateParameters:
            print self.rateParameters
            outfile.write('\n')

            nCR= len(self.CR)
            print "length of CR: ", nCR
            #print self.regionMapping
            shift_SR=0

            if shapeFile :
                for i_cr, cr in enumerate(self.CR):
                        if "And" in p[0] :
                            outfile.write('%s_norm_%s rateParam %s %s %s %s\n'%("Bin{}".format(i_cr), p[0], "Bin{}".format(i_cr), p[0].split("And")[0], str(p[1]), str(p[2]) ))
                            outfile.write('%s_norm_%s rateParam %s %s %s %s\n'%("Bin{}".format(i_cr), p[0], "Bin{}".format(i_cr), p[0].split("And")[1], str(p[1]), str(p[2]) ))
                        else :
                            outfile.write('%s_norm_%s rateParam %s %s %s %s\n'%("Bin{}".format(i_cr), p[0], "Bin{}".format(i_cr), p[0], str(p[1]), str(p[2]) ))
                        
                        
            else :
                for i_cr, cr in enumerate(self.CR):
                        outfile.write('%s_norm_%s rateParam %s %s (@0*1) %s_norm_%s\n'%("Bin{}".format(i_cr), p[0], "Bin{}".format(i_cr), p[0],"Bin{}".format(i_cr), p[0] ))
                        print '%s_norm_%s rateParam %s %s (@0*1) %s_norm_%s\n'%("Bin{}".format(i_cr), p[0], "Bin{}".format(i_cr), p[0],"Bin{}".format(i_cr), p[0] )
                        corr_SR = self.SR[shift_SR:shift_SR+self.regionMapping[i_cr]]
                        
                        for i_sr, sr in enumerate(corr_SR):
                            print '%s_norm_%s rateParam %s %s (@0*1) %s_norm_%s\n'%("Bin{}".format(nCR+shift_SR+i_sr), p[0], "Bin{}".format(nCR+shift_SR+i_sr), p[0],"Bin{}".format(i_cr), p[0] )
                            outfile.write('%s_norm_%s rateParam %s %s (@0*1) %s_norm_%s\n'%("Bin{}".format(nCR+shift_SR+i_sr), p[0], "Bin{}".format(nCR+shift_SR+i_sr), p[0],"Bin{}".format(i_cr), p[0] ))
                        
                        shift_SR +=len(corr_SR)
                        # outfile.write('%s_norm_%s extArg %s %s\n'%("Bin{}".format(i_cr), p[0], str(p[1]), str(p[2])))
		
        for p in self.freeParameters:
            outfile.write('\n')
            #for b in self.freeParamBins[p[0]]:
	    #bins = [val for keys,val in region.items() if p[0] in names]

            for b in self.bins:
                outfile.write('%s rateParam %s %s %s %s\n'%(p[0], b, p[1], str(p[2]), str(p[3])))
            #for b in self.bins:
            #    outfile.write('%s rateParam %s %s %s %s\n'%(p[0], b, p[1], str(p[2]), str(p[3])))

        if shapeFile and not noMCStat:
            outfile.write('* autoMCStats 10 \n')

        outfile.close()
        print "[cardFileWrite] Written card file %s"%fname
        return fname

    def makeHist(self, name):
        return ROOT.TH1F(name, name, len(self.bins), 0, len(self.bins))


    def getPTgroupedHistos(self, hist_object) :
        nCR= len(self.CR)
        h_pT_grouped = []
        shift_SR=0
        for i_cr, cr in enumerate(self.CR):
            nbins = self.regionMapping[i_cr] + 1
            
            firstBin = 0 
            h_pT_grouped.append(ROOT.TH1F("Bin{}_{}".format(i_cr,hist_object.GetName()), "Bin{}_{}".format(i_cr,hist_object.GetName()), nbins, 0, nbins))
            # h_pT_grouped.append(ROOT.TH1F("Bin{}_{}".format(i_cr,hist_object.GetName()), "Bin{}_{}".format(i_cr,hist_object.GetName()), binning[nbins]))
            

            corr_SR = self.SR[shift_SR:shift_SR+self.regionMapping[i_cr]]
           
            for i_sr, sr in enumerate(corr_SR) :
                h_pT_grouped[i_cr].SetBinContent(firstBin+1,hist_object.GetBinContent(nCR+shift_SR+i_sr+1)) #12+0+{0,1,2,3}
                h_pT_grouped[i_cr].SetBinError(firstBin+1,hist_object.GetBinError(nCR+shift_SR+i_sr+1)) 
                
                firstBin += 1
            shift_SR += len(corr_SR)


            h_pT_grouped[i_cr].SetBinContent(firstBin+1,hist_object.GetBinContent(i_cr+1))
            h_pT_grouped[i_cr].SetBinError(firstBin+1,hist_object.GetBinError(i_cr+1))
            
        return h_pT_grouped


    def writeToShapeFile(self, fname, noMCStat=False):
        bins        = natural_sort(self.bins)
        processes   = []
        nuisances   = [ u for u in self.uncertainties if (not 'stat' in u.lower() and self.uncertaintyString[u] == 'shape')  ] # stat uncertainties are treated differently
        logNormal   = [ u for u in self.uncertainties if (not 'stat' in u.lower() and self.uncertaintyString[u] == 'lnN')  ]
        
        
        for b in bins:
            for p in self.processes[b]:
                if p not in processes: processes.append(p)

        # define a dict that stores which shape nuisances are relevant for each process
        nuisForProc = {proc:[] for proc in processes}

        # first, fill the observations. easy.
        data_obs = self.makeHist('data_obs')
        for i,b in enumerate(bins):
            data_obs.SetBinContent(i+1, self.observation[b])

        # predictions, using the correct stats
        histos = {}
        for process in processes:
            histos[process] = self.makeHist(process)
            
            # create the histograms for all uncertainties that are relevant for the process
            for unc in nuisances:
                for i,b in enumerate(bins):
                    if self.uncertaintyVal.has_key((unc, b, process)):
                        if not (unc.lower().count('up') or unc.lower().count('down')):
                            nuisForProc[process].append(unc)
                            up      = '%s_%sUp'%(process, unc)
                            down    = '%s_%sDown'%(process, unc)
                            try:
                                histos[up]
                            except:
                                histos[up]   = self.makeHist(up)
                                histos[down] = self.makeHist(down)
                        else:
                            if not unc in nuisForProc[process]:
                                nuisForProc[process].append(unc.replace('Up','').replace('Down',''))
                            try:
                                histos['%s_%s'%(process, unc)]
                            except:
                                histos['%s_%s'%(process, unc)] = self.makeHist('%s_%s'%(process, unc))

            # fill the histograms of central values and uncertainties. if no uncertainty exists for a bin the up/down variations are set to the central values as well.
            for i,b in enumerate(bins):
                if self.expectation.has_key((b, process)):
                    expect  = self.expectation[(b, process)]
                    try:
                        relUnc  = self.uncertaintyVal[('Stat_'+b+'_'+process, b, process)]
                    except:
                        relUnc  = 1
                    unc     = (relUnc-1)*expect
                    histos[process].SetBinContent(i+1, expect)
                    histos[process].SetBinError(i+1, unc)
                    for unc in nuisances:
                        if self.uncertaintyVal.has_key((unc, b, process)):
                            relUnc  = self.uncertaintyVal[(unc, b, process)]
                        else:
                            relUnc = 1.
                        if unc in nuisForProc[process]:
                            if not (unc.lower().count('up') or unc.lower().count('down')):
                                histos['%s_%sUp'%(process, unc)].SetBinContent(i+1, relUnc*expect)
                                # histos['%s_%sDown'%(process, unc)].SetBinContent(i+1, (2-relUnc)*expect)
                                histos['%s_%sDown'%(process, unc)].SetBinContent(i+1, expect/relUnc)
                            else:
                                histos['%s_%s'%(process, unc)].SetBinContent(i+1, relUnc*expect)
        # define the file names
        rootFile = fname.split('/')[-1]
        rootFileFull = fname
        txtFile  = fname.replace('.root', 'Card.txt')

        # create the card-file with and group all CR and SR bins together
        self.bins = self.bins[0:self.nCR]

        # self.niceNames[bins[0]] = 'inclusive bin'
        shift_SR = 0
        for i_bin in xrange(self.nCR) :
            s_start = self.nCR + shift_SR 
                
            for process in processes:
                # s_start = self.nCR + shift_SR 
                self.specifyExpectation(self.bins[i_bin], process, histos[process].GetBinContent(i_bin+1) + histos[process].Integral(s_start+1,s_start+self.regionMapping[i_bin]))
                
                

                #     exit(0)
                #print "setting expectation to shape file", self.bins[0], process,  histos[process].Integral()
                #need to set the strength of the shape uncertainties defined by the histograms. If not turned on (meaning a value greater than 0), they have no effect.
                for nuis in nuisForProc[process]:
                    self.specifyUncertainty(nuis, bins[i_bin], process, 1)
            
            shift_SR += self.regionMapping[i_bin]
                
            self.specifyObservation(self.bins[i_bin], int(round(data_obs.GetBinContent(i_bin+1) + data_obs.Integral(s_start+1,s_start+self.regionMapping[i_bin]),0) )  )

        self.uncertainties = logNormal + nuisances # need to fix things if up/down are already provided
        self.writeToFile(txtFile, shapeFile=rootFile, noMCStat=noMCStat)
        







        # write all the histograms to a root file
        h_pT_grouped = self.getPTgroupedHistos(data_obs)
        for ihh, hhhisto in enumerate(h_pT_grouped) :
            if ihh == 0 :
                writeObjToFile(rootFileFull, hhhisto)
            else :
                writeObjToFile(rootFileFull, hhhisto, update=True)

        for h in sorted(histos.keys()):
            h_pT_grouped = self.getPTgroupedHistos(histos[h])
            for ihh, hhhisto in enumerate(h_pT_grouped) :
                writeObjToFile(rootFileFull, hhhisto, update=True)
        
        
        return txtFile

    def combineCards(self, cards):

        import uuid, os
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join(self.releaseLocation, ustr)
        logger.info("Creating %s", uniqueDirname)
        os.makedirs(uniqueDirname)

        years = cards.keys()
        cmd = ''
        for year in years:
            cmd += " dc_%s=%s"%(year, cards[year])

        combineCommand  = "cd "+uniqueDirname+";combineCards.py %s > combinedCard.txt; text2workspace.py combinedCard.txt --X-allow-no-signal -m 125"%(cmd)
        os.system(combineCommand)
        resFile = cards[years[0]].replace(str(years[0]), 'COMBINED')
        f = resFile.split('/')[-1]
        resPath = resFile.replace(f, '')
        if not os.path.isdir(resPath):
            os.makedirs(resPath)
        logger.info("Putting combined card into dir %s", resPath)
        shutil.copyfile(uniqueDirname+"/combinedCard.txt", resFile)
        shutil.copyfile(uniqueDirname+"/combinedCard.root", resFile.replace(".txt",".root"))

        return resFile

    def readResFile(self, fname):
        import ROOT
        f = ROOT.TFile.Open(fname)
        t = f.Get("limit")
        l = t.GetLeaf("limit")
        qE = t.GetLeaf("quantileExpected")
        limit = {}
        preFac = 1.
        for i in range(t.GetEntries()):
                t.GetEntry(i)
                limit["{0:.3f}".format(round(qE.GetValue(),3))] = preFac*l.GetValue()
        f.Close()
        return limit

    def calcLimit(self, fname=None, options="", normList=[]):
        import uuid, os
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join(self.releaseLocation, ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)

        if fname is not None:  # Assume card is already written when fname is not none
          filename = os.path.abspath(fname)
        else:
          filename = fname if fname else os.path.join(uniqueDirname, ustr+".txt")
          self.writeToFile(filename)
        resultFilename = filename.replace('.txt','')+'.root'

        assert os.path.exists(filename), "File not found: %s"%filename
        
        if normList:
            from Analysis.Tools.cardFileWriter.getNorms import getNorms
            filename += " > output.txt"
        
        combineCommand = "cd "+uniqueDirname+";combine --saveWorkspace -M AsymptoticLimits %s %s"%(options,filename)
        print combineCommand
        os.system(combineCommand)

        tempResFile = uniqueDirname+"/higgsCombineTest.AsymptoticLimits.mH120.root"

        try:
            res= self.readResFile(tempResFile)
        except:
            res=None
            print "[cardFileWrite] Did not succeed reeding result."
        if res:
            shutil.copyfile(tempResFile, resultFilename)
        
        if normList:
            getNorms(dirName=uniqueDirname, normsToExtract=normList)
            shutil.copyfile(uniqueDirname + '/SF.txt', resultFilename.replace('.root','_SF.txt'))
        
        shutil.rmtree(uniqueDirname)
        return res

    def readNLLFile(self, fname):
        import ROOT
        f = ROOT.TFile.Open(fname)
        t = f.Get("limit")
        nll = {}
        t.GetEntry(0) # changed from 1!
        # prefit NLL
        nll["nll0"] = t.nll0
        # delta NLL to prefit (should always be negative since stuff is fitted)
        nll["nll"] = t.nll
        # absolute NLL postfit
        nll["nll_abs"] = t.nll0 + t.nll
        f.Close()
        return nll

    def calcNLL(self, fname=None, options=""):
        '''
        Does max likelihood fits, both with r=1 and a best-fit value
        '''
        import uuid, os
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join(self.releaseLocation, ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)
        if fname is not None:  # Assume card is already written when fname is not none
          filename = os.path.abspath(fname)
        else:
          filename = fname if fname else os.path.join(uniqueDirname, ustr+".txt")
          self.writeToFile(filename)

        combineCommand  = "cd "+uniqueDirname+";eval `scramv1 runtime -sh`;combine -M MultiDimFit -n Nominal --saveNLL --forceRecreateNLL --freezeParameters r %s %s"%(options,filename)
        os.system(combineCommand)
        nll = self.readNLLFile(uniqueDirname+"/higgsCombineNominal.MultiDimFit.mH120.root")
        nll["bestfit"] = nll["nll"]
        shutil.rmtree(uniqueDirname)

        return nll

    def poissonPrefitNLL( self ):
        '''
        Compute the prefit NLL from expectation and observation without using combine
        '''
        obs_sorted = [self.observation[b] for b in self.bins]
        exp_sorted = map( lambda l: sum([self.expectation[k] for k in l]), [filter( lambda k:k[0]==bin, self.expectation.keys() ) for bin in self.bins ])
        return -sum([ -lam + obs*log(lam) - sum( [log(n_) for n_ in range(1, obs+1 )] ) for obs, lam in zip(obs_sorted, exp_sorted)])

    def calcNuisances(self, fname=None, options="", outputFileAddon = "", bonly=False):
        import uuid, os
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join(self.releaseLocation, ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)
        shutil.copyfile(os.path.join(os.environ['CMSSW_BASE'], 'src', 'Analysis', 'Tools', 'python', 'cardFileWriter', 'diffNuisances.py'), os.path.join(uniqueDirname, 'diffNuisances.py'))

        if fname is not None:  # Assume card is already written when fname is not none
          filename = os.path.abspath(fname)
        else:
          filename = fname if fname else os.path.join(uniqueDirname, ustr+".txt")
          self.writeToFile(filename)
        resultFilename      = filename.replace('.txt','')+'_nuisances.txt'
        resultFilenameFull  = filename.replace('.txt','')+'_nuisances_full.txt'
        resultFilename2     = filename.replace('.txt','')+'_nuisances.tex'
        resultFilename2Full = filename.replace('.txt','')+'_nuisances_full.tex'

        assert os.path.exists(filename), "File not found: %s"%filename

        combineCommand  = "cd "+uniqueDirname+";combine --robustHesse 1 --forceRecreateNLL -M FitDiagnostics --saveShapes --saveNormalizations --saveOverall --saveWithUncertainties %s %s"%(options,filename)
        combineCommand +=";python diffNuisances.py  fitDiagnosticsTest.root &> nuisances.txt"
        combineCommand +=";python diffNuisances.py -a fitDiagnosticsTest.root &> nuisances_full.txt"
        if bonly:
          combineCommand +=";python diffNuisances.py -bf latex fitDiagnosticsTest.root &> nuisances.tex"
          combineCommand +=";python diffNuisances.py -baf latex fitDiagnosticsTest.root &> nuisances_full.tex"
        else:
          combineCommand +=";python diffNuisances.py -f latex fitDiagnosticsTest.root &> nuisances.tex"
          combineCommand +=";python diffNuisances.py -af latex fitDiagnosticsTest.root &> nuisances_full.tex"
        print combineCommand
        os.system(combineCommand)

        if outputFileAddon: outputFileAddon = "_"+outputFileAddon
        shutil.copyfile(uniqueDirname+'/fitDiagnosticsTest.root', fname.replace('.txt','%s_FD.root'%(outputFileAddon)))

        tempResFile      = uniqueDirname+"/nuisances%s.txt"%(outputFileAddon)
        tempResFileFull  = uniqueDirname+"/nuisances%s_full.txt"%(outputFileAddon)
        tempResFile2     = uniqueDirname+"/nuisances%s.tex"%(outputFileAddon)
        tempResFile2Full = uniqueDirname+"/nuisances%s_full.tex"%(outputFileAddon)
        shutil.copyfile(tempResFile, resultFilename)
        shutil.copyfile(tempResFileFull, resultFilenameFull)
        shutil.copyfile(tempResFile2, resultFilename2)
        shutil.copyfile(tempResFile2Full, resultFilename2Full)

        shutil.rmtree(uniqueDirname)
        return


    def calcSignif(self, fname="", options=""):
        import uuid, os
        ustr          = str(uuid.uuid4())
        uniqueDirname = os.path.join(self.releaseLocation, ustr)
        print "Creating %s"%uniqueDirname
        os.makedirs(uniqueDirname)
        combineCommand = "cd "+uniqueDirname+";combine --saveWorkspace -M ProfileLikelihood --uncapped 1 --significance --rMin -5 %s %s"%(options,fname)
        os.system(combineCommand)
        try:
            res= self.readResFile(uniqueDirname+"/higgsCombineTest.ProfileLikelihood.mH120.root")
            os.system("rm -rf "+uniqueDirname+"/higgsCombineTest.ProfileLikelihood.mH120.root")
        except:
            res=None
            print "Did not succeed."
        shutil.rmtree(uniqueDirname)

        return res
