'''
Some definitions.
'''
import os
if os.environ['USER'] in ['llechner']:
    plot_directory         = "/afs/hephy.at/user/l/llechner/www/TTGammaEFT/"
    cache_directory        = "/afs/hephy.at/data/llechner01/TTGammaEFT/cache/"
    cern_proxy_certificate = "/afs/cern.ch/user/l/llechner/private/.proxy"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/llechner/'
elif os.environ['USER'] in ['lukas.lechner']:
    plot_directory         = "/mnt/hephy/cms/lukas.lechner/www/TTGammaEFT/"
    cache_directory        = "/users/lukas.lechner/public/cache/"
#    cache_directory        = "/mnt/hephy/cms/lukas.lechner/TTGammaEFT/cache/"
    cern_proxy_certificate = "/users/lukas.lechner/private/.proxy"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/llechner/'
elif os.environ['USER'] in ['rosmarie.schoefbeck']:
    plot_directory         = "/mnt/hephy/cms/rosmarie.schoefbeck/www/TTGammaEFT/"
    cache_directory        = "/users/rosmarie.schoefbeck/public/cache/"
    cern_proxy_certificate = "/users/rosmarie.schoefbeck/private/.proxy"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/llechner/'
elif os.environ['USER'] in ['janik.andrejkovic']:
    plot_directory         = "/mnt/hephy/cms/janik.andrejkovic/www/"
    cache_directory        = "/users/janik.andrejkovic/public/tWZ/cache/"
    cern_proxy_certificate = "/users/janik.andrejkovic/private/.proxy"
elif os.environ['USER'] in ['dspitzba', 'dspitzbart']:
    plot_directory         = "/afs/hephy.at/user/d/dspitzbart/www/stopsDileptonLegacy/"
    cache_directory        = "/afs/hephy.at/data/dspitzbart01/cache/"
    cern_proxy_certificate = "/afs/cern.ch/user/d/dspitzba/private/.proxy"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/dspitzbart/'
elif os.environ['USER'] in ['phussain']:
    cache_directory        = "/afs/hephy.at/data/dspitzbart01/cache/"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/prhussai/'
elif os.environ['USER'] in ['priya.hussain']:
    cache_directory        = "/mnt/hephy/cms/priya.hussain/StopsCompressed/cache/"
    plot_directory         = "/mnt/hephy/cms/priya.hussain/www/StopsCompressed/"
elif os.environ['USER'] in ['dietrich.liko']:
    cache_directory        = "/groups/hephy/cms/dietrich.liko/StopsCompressed/cache/"
    plot_directory         = "/groups/hephy/cms/dietrich.liko/www/StopsCompressed/"
elif os.environ['USER'] in ['rschoefbeck']:
    plot_directory         = "/afs/hephy.at/user/r/rschoefbeck/www/StopsDilepton/"
    cache_directory        = "/afs/hephy.at/data/rschoefbeck01/cache/"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/schoef/'
elif os.environ['USER'] in ['robert.schoefbeck']:
    plot_directory         = "/groups/hephy/cms/robert.schoefbeck/www/Analysis/"
    cache_directory        = "/groups/hephy/cms/robert.schoefbeck/caches/"
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/schoef/'
    remote_host            = 'schoef@lxplus.cern.ch'
    remote_www_directory   = '/eos/user/s/schoef/www'
elif os.environ['USER'] in ['lena.wild']:
    plot_directory         = "/groups/hephy/cms/lena.wild/www/Analysis/"
    cache_directory        = "/groups/hephy/cms/lena.wild/caches/"
    remote_host            = 'lwild@lxplus.cern.ch'
    remote_www_directory   = '/eos/user/l/lwild/www'
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
elif os.environ['USER'] in ['dennis.schwarz']:
    plot_directory         = "/mnt/hephy/cms/dennis.schwarz/www/tWZ/"
    cache_directory        = "/users/dennis.schwarz/public/cache/"
    cern_proxy_certificate = "/users/dennis.schwarz/private/.proxy"
elif os.environ['USER'] in ['cristina.giordano']:
    plot_directory         = "/groups/hephy/cms/cristina.giordano/www/tttt/"
    cache_directory        = "/users/cristina.giordano/public/cache/"
    cern_proxy_certificate = "/users/cristina.giordano/private/.proxy"
else:
    plot_directory         = "/afs/hephy.at/user/%s/%s/www/Analysis/"%(os.environ['USER'][0],os.environ['USER'])
    cache_directory        = "/afs/hephy.at/data/%s01/Analysis/"%os.environ['USER']
    cern_proxy_certificate = "/afs/cern.ch/user/%s/%s/private/.proxy"%(os.environ['USER'][0],os.environ['USER'])
    dpm_directory          = '/dpm/oeaw.ac.at/home/cms/store/user/%s/'%os.environ['USER']

if not os.path.isdir(cache_directory): os.makedirs(cache_directory)
