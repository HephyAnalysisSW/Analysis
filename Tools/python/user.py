'''
Some definitions.
'''
import os
if os.environ['USER'] in ['llechner']:
    plot_directory         = "/afs/hephy.at/user/l/llechner/www/TTGammaEFT/"
    cache_directory        = "/afs/hephy.at/data/llechner01/TTGammaEFT/cache/"
    cern_proxy_certificate = "/afs/cern.ch/user/l/llechner/private/.proxy"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/llechner/'
elif os.environ['USER'] in ['dspitzba', 'dspitzbart']:
    plot_directory         = "/afs/hephy.at/user/d/dspitzbart/www/stopsDileptonLegacy/"
    cache_directory        = "/afs/hephy.at/data/dspitzbart01/cache/"
    cern_proxy_certificate = "/afs/cern.ch/user/d/dspitzba/private/.proxy"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/dspitzbart/'
elif os.environ['USER'] in ['phussain']:
    cache_directory        = "/afs/hephy.at/data/dspitzbart01/cache/"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/prhussai/'
elif os.environ['USER'] in ['rschoefbeck']:
    plot_directory         = "/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/"
    cache_directory        = "/afs/hephy.at/data/rschoefbeck01/cache/"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/schoef/'
elif os.environ['USER'] in ['schoef']:
    plot_directory         = "/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/"
    cache_directory        = "/afs/hephy.at/data/rschoefbeck01/cache/"
    cern_proxy_certificate = "/afs/cern.ch/user/s/schoef/private/.proxy"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/schoef/'
elif os.environ['USER'] in ['ttschida']:
    plot_directory         = "/afs/hephy.at/user/t/ttschida/www/StopsDilepton/"
    cache_directory        = "/afs/hephy.at/data/cms04/ttschida/cache/"
    cern_proxy_certificate = "/afs/cern.ch/user/t/ttschida/private/.proxy"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/ttschida/'
else:
    plot_directory         = "/afs/hephy.at/user/%s/%s/www/Analysis/"%(os.environ['USER'][0],os.environ['USER'])
    cache_directory        = "/afs/hephy.at/data/%s01/Analysis/"%os.environ['USER']
    cern_proxy_certificate = "/afs/cern.ch/user/%s/%s/private/.proxy"%(os.environ['USER'][0],os.environ['USER'])
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/%s/'%os.environ['USER']

if not os.path.isdir(cache_directory): os.makedirs(cache_directory)
