#!/usr/bin/env python
import os
from Analysis.Tools.helpers import  checkRootFile
from subprocess import call

def get_parser():
    ''' Argument parser for post-processing module.
    '''
    import argparse 
    argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")

    argParser.add_argument('--dir',
        action='store',
        nargs='?',
        type=str,
        default='/scratch/rschoefbeck/test',
        help="Name of the directory to be re-hadded"
        )

    argParser.add_argument('--skip',
        action='store',
        nargs='*',
        type=str,
        default=[ 'T2tt', 'TTbarDMJets' ],
        help="Substrings in filename that should not be touched."
        )

#    argParser.add_argument('--targetDir',
#        action='store',
#        nargs='?',
#        type=str,
#        default=None,
#        help="Where to store output. Default is in the same place."
#        )

    argParser.add_argument('--treeName',
        action='store',
        nargs='?',
        type=str,
        default='Events',
        help="Treename."
        )

    argParser.add_argument('--sizeGB',
        action='store',
        nargs='?',
        type=float,
        default=2.,
        help="What size you want?"
        )

    argParser.add_argument('--really',
        action='store_true',
        help="Just to be sure.")

    argParser.add_argument('--delete',
        action='store_true',
        help="Delete originals?")

    argParser.add_argument('--overwrite',
        action='store_true',
        help="Delete originals?")

    argParser.add_argument('--logLevel', 
        action='store',
        nargs='?',
        choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],
        default='INFO',
        help="Log level for logging"
        )
        
    return argParser

options = get_parser().parse_args()

if options.treeName=='': 
    options.treeName=None

# Logger
import Analysis.Tools.logger as logger
logger = logger.get_logger(options.logLevel, logFile = None )

# Walk the directory structure and group files in 'jobs' of [f1_0.root, f1_1.root, ...]  tootalling to approx. sizeGB
jobs = []
for dirName, subdirList, fileList in os.walk(options.dir):
    rootFiles = []
    for f in fileList:
        if f.endswith('.root'):
            full_filename = os.path.join(dirName, f)
            if not '_reHadd_' in f:
                to_skip = False
                for skip in options.skip:
                    if skip in f:
                        logger.info( "Found skip string %s in %s. Skipping.", skip, f )
                        to_skip = True
                        break
                if to_skip: continue
                isOK =  checkRootFile( full_filename, checkForObjects = [options.treeName]) \
                        if options.treeName is not None else checkRootFile( full_filename )
                if isOK:
                    rootFiles.append( f )
                else:
                    logger.warning( "File %s does not look OK. Checked for tree: %r", full_filename, options.treeName )
            else:
                logger.info( "Found '_reHadd_' in file %s in %s. Skipping.", full_filename, dirName )
    job = []
    jobsize = 0
    for fname in rootFiles:
        filename, file_extension = os.path.splitext(fname)
        n_str = filename.split('_')[-1]
        if n_str.isdigit():
            full_filename = os.path.join(dirName, fname)
            jobsize += os.path.getsize( full_filename  )
            job.append( full_filename )
            if jobsize>1024**3*options.sizeGB:
                jobs.append(job)
                job = []
                jobsize = 0
        else:
            logger.warning( "Could not determine file counter for file %s in %s. Expected: fileName_<int>.root", fname, dirName )

    if len(job)>0: jobs.append(job)

basename_counter = {}
for job in jobs:

    fnames = list(set([ '_'.join( os.path.splitext(fname)[0].split('_')[:-1] ) for fname in job]))
    if not len(fnames)==1:
        logger.warning("Problem in job %r: Filenames are not consistent. Found %r. Do nothing.", job, fnames)
        continue
    basename = fnames[0]
    if not basename_counter.has_key(basename): 
        basename_counter[basename] = 0
    else:
        basename_counter[basename] += 1

    targetFileName = basename+'_reHadd_'+str(basename_counter[basename])+'.root'
    size = sum(os.path.getsize(f) for f in job)/1024.**3
    logger.info( "Hadding %i files to %s to a total of %3.2f GB.", len(job), targetFileName, size )

    if targetFileName in job:
        raise ValueError( "Found file %s in job %r. Should not happen and can't hadd this way.", targetFileName, job )

    if options.really:
        cmd = ["hadd"]
    else:
        cmd = ["echo", "hadd"]
#        logger.info("Use --really to do it really. Use --delete to delete the input files." )

    if os.path.exists( targetFileName ):
        if options.overwrite:
            logger.info( "Found file %s. Overwriting.", targetFileName )
            cmd.append( '-f' )
        else:
            logger.info( "File %s already exists. Skipping.", targetFileName )
            cmd = None
   
    from subprocess import call

    if cmd is not None: call(cmd + [targetFileName] + job) 
    
    if options.delete:
        isOK =  checkRootFile( targetFileName, checkForObjects = [options.treeName]) \
            if options.treeName is not None else checkRootFile( targetFileName ) 
        if isOK:
            for f in job:
                os.remove( f )
                logger.debug( "Deleted input %s", f)
        else:
            logger.warning( "File %s does not look OK. Checked for tree: %r. Did not delete input.", targetFileName, options.treeName )
    else:
        if not isOK: logger.warning( "File %s does not look OK. Checked for tree: %r.", targetFileName, options.treeName )
            

    if options.really:
        logger.info("Done.")
    else:
        logger.info("Done with nothing. Use --really to hadd and --delete to delete the input.")
 
