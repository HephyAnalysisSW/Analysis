'''
find datasets for files on EOS user directory
Run 
find /eos/vbc/experiments/cms/store/user/ -name "*.root" >& all_files.txt &
and 
ipython -i find_samples.py -- all_files.txt
'''

import sys
import os
import re

patterns = {}

with open(sys.argv[1]) as fp:
    lines = fp.readlines()
    for line in lines:
        #if not 'schoef' in line:continue
        line = line.strip()
        f_name = re.findall( '/[0-9]+/[A-Za-z0-9\-_]*\.root', line)
        if not f_name:
            continue
        f_name = f_name[0]
        pattern = line.replace(f_name,'')
        # drop last three folder levels in case of data
        if "/store/data/" in pattern:
            s = pattern.rstrip("/").split("/")
            pattern = "/".join(s[:s.index('store')+6])
        patterns[pattern] = '/store/'+line.split('/store/')[-1]

def get_size(start_path = '.'):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            # skip if it is symbolic link
            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)

    return total_size

def read_from_das( job ):
    pattern, f = job
    
    #dbs='dasgoclient -query="dataset file=%s instance=prod/phys03"'%f
    dbs='dasgoclient -query="dataset file=%s"'%f
    print "Query:", dbs
    dbsOut = os.popen(dbs)
    lines = dbsOut.readlines()

    files = []
    result = {'name':"not_published", 'size':get_size(pattern), 'location': pattern}
    for line in lines:
        if line.startswith('/'):
            #line = line.rstrip()
            #filename = redirector+'/'+line
            result['name'] = line.rstrip()
            print "Found", result 
            break
    dbsOut.close()
    return result 

from multiprocessing import Pool

jobs = list(patterns.iteritems())
p = Pool(30)
results = p.map(read_from_das, jobs)

consumption = {}
for dict_ in results:
    user = 'unknown'
    try:
        user = dict_['location'].split('/cms/store/user/')[-1].split('/')[0]
    except:
        pass
    if not consumption.has_key(user):
        consumption[user] = {'size':0, 'samples' : []}
    consumption[user]['size']+=dict_['size']
    consumption[user]['samples'].append(dict_)

def human_size(bytes, units=[' bytes','kB','MB','GB','TB', 'PB', 'EB']):
    """ Returns a human readable string representation of bytes """
    return str(bytes) + units[0] if bytes < 1024 else human_size(bytes>>10, units[1:])
 
u_list = consumption.keys()
u_list.sort(key=lambda u:-consumption[u]['size']) 
with open( 'datasets.txt', 'w') as f:
    for user in u_list:
        f.write("%8s total consumption: %s\n"%(user, human_size(consumption[user]['size'])))
    for user in u_list:
        f.write('\n')
        f.write("%s total consumption: %s\n"%(user, human_size(consumption[user]['size'])))
        consumption[user]['samples'].sort(key=lambda s:-s['size'])
        for sample in consumption[user]['samples']:
            f.write("  %10s DAS: %s EOS: %s\n"%( human_size(sample['size']), sample['name'], sample['location'])) 
        f.write('\n')

#    for directory, dataset in sorted(datasets.iteritems()):
#        #print directory, dataset
#        f.write('%s %s\n'% (directory, dataset))

