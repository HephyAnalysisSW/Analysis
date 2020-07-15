#Also works the other way around! 
import Analysis.Tools.syncer 
import os

import Analysis.Tools.logger as logger
logger  = logger.get_logger('DEBUG', logFile = None)

if not os.path.isdir('www'):
    os.makedirs( 'www')

# ROOT example
import ROOT
c1 = ROOT.TCanvas()
c1.Print('www/x.png')
c1.Print('www/y.pdf')

#pickle example
import pickle
x = {}
pickle.dump( 'x', file('www/z.pkl','w'))
