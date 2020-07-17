#!/usr/bin/env python
import os, sys
import subprocess

redir       = "root://hephyse.oeaw.ac.at/"
hostname    = os.getenv("HOSTNAME")


def checkRootFile( file ):
    logger.info("Checking root file: %s"%file)
    from Analysis.Tools.helpers import checkRootFile, deepCheckRootFile, deepCheckWeight
    valid = checkRootFile( file, checkForObjects=["Events"] ) and deepCheckRootFile( file ) and deepCheckWeight( file )
    if valid:
        logger.info("Check done!")
    else:
        logger.info("Corrupt root file: %s"%file)
    return valid

def getDPMFiles( path, fromLxPlus=False ):
    """ get all files in dpm directory
    """
    if "cern" in hostname or "clip" in hostname:
        cmd = "xrdfs %s ls %s" %(redir,path)
    else:
        cmd = "dpns-ls %s"%path
    p = subprocess.Popen( [cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    return [ item.rstrip("\n") for item in p.stdout.readlines() ]

def convertToPathList( paths, local=False ):
    """ convert path with ending * to list of pathes with all possible files
    """
    if not isinstance(paths, list) and not isinstance(paths, basestring):
        print "convertToPathList needs a string or list as input. Exiting..."
        sys.exit(0)

    if isinstance(paths, basestring) and not "*" in paths: return [paths]
    elif isinstance(paths, basestring): paths = [paths]

    allPathList = []
    for path in paths:
        pathList = []
#        print "Checking path " + path
        dirs = path.split("/")
        dirPath = ""
        for i, dir in enumerate(dirs):
            if not "*" in dir:
                dirPath += dir + "/"
                continue
#            print "Checking dirPath " + dirPath
            allFiles = os.listdir(dirPath) if local else getDPMFiles(dirPath)
            for file in allFiles:
#                print "Checking file " + file
                subString = file
                copyFile  = True
                # check if all parts are contained in filename in the right order
                for part in dir.split("*"):
                    if not part: continue
                    if not part in subString:
                        copyFile = False
                        break
                    subString = "".join( subString.split( part )[1:] )
                if copyFile:
                    file = "/".join( [file]+dirs[i+1:] )
#                    pathList.append( os.path.join( dirPath, file ) )
#                    print
#                    print "Found path to print: " + os.path.join( dirPath, file )
#                    print
                    pathList += convertToPathList( os.path.join( dirPath, file ), local=local )
#            if copyFile: break       
            if pathList: break       
#        print pathList
        if not pathList and not "*" in path: allPathList += [path]
        elif not pathList and "*" in path: continue
        else: allPathList += pathList
    return allPathList

def removeDPMFiles( path ):
    """ remove dpm dir or file (only for user dirs)
    """

    if not os.environ["USER"] in path and not (os.environ["CERN_USER"] in path and "clip" in hostname):
        logger.info( "ATTENTION: DO NOT REMOVE FILES FROM OTHERS DPM DIRECTORIES: %s"%path )
        logger.info( "EXITING" )
        sys.exit(1)

    rmPathList = convertToPathList( path )

    for rmPath in rmPathList:
        logger.info( "Removing %s"%rmPath )
        if "cern" in hostname or "clip" in hostname:
            try:    os.system( "xrdfs %s rmdir %s" %(redir,rmPath) )
            except: os.system( "xrdfs %s rm %s" %(redir,rmPath) )
        else:
            os.system( "/usr/bin/rfrm -rf %s"%rmPath )

def checkDPMDirExists( path ):

    if path.endswith("*"):
        logger.info( "Cannot create directory: %s"%path )
        sys.exit(1)

    checkDir  = getDPMFiles( path )
    dirExists = len( [ item for item in checkDir if not "No such file or directory" in item and not "error" in item ] ) > 0
    return dirExists

def makeDPMDir( path ):
    """ create a dpm dir
    """

    if checkDPMDirExists( path ):
        logger.info( "Directory exists: %s"%path )
    else:
        motherPath = "/".join( path.split("/")[:-1] )
        if not checkDPMDirExists( motherPath ):
            makeDPMDir( motherPath )
        logger.info( "Creating directory: %s"%path )
        if "cern" in hostname or "clip" in hostname:
            cmd = "xrdfs %s mkdir %s" %(redir,path)
        else:
            cmd = "dpns-mkdir %s"%path
        os.system( cmd )        

def isFile( path ):
    """ check if path is a file
    """

    return "." in path.split("/")[-1]

def copyDPMFiles( fromPath, toPath, toLocal=False, fromLocal=False, rootFileCheck=False ):
    """ copy files or directories including subdirectories
    """

    fromPathList = convertToPathList( fromPath, local=fromLocal )
    if toLocal:
        if not os.path.isdir( toPath ): os.makedirs( toPath )
    else:
        makeDPMDir( toPath )

    for file in fromPathList:
        if isFile( file ):
            logger.info( "Copying %s to %s"%( file, toPath ) )
            cmd = "xrdcp -r %s%s %s%s"%(redir if not fromLocal else "", file, redir if not toLocal else "", toPath)
            os.system( cmd )

            if rootFileCheck:
                corrupt = True
                for i in range(10):
                    if checkRootFile( "%s%s/%s"%(redir if not toLocal else "", toPath, file.split("/")[-1]) ):
                        corrupt = False
                        break
                    # remove files if corrupt (dont do -rf, otherwise it always copies)
                    logger.info("Corrupt root file! Removing file and trying again!")
                    if toLocal: os.system("rm %s/%s"%(toPath,file.split("/")[-1]))
                    else:       removeDPMFiles( "%s/%s"%(toPath,file.split("/")[-1]) )
                    os.system( cmd ) #try again if corrupt
                if corrupt:
                    if toLocal: os.system("rm %s/%s"%(toPath,file.split("/")[-1]))
                    else:       removeDPMFiles( "%s/%s"%(toPath,file.split("/")[-1]) )

        else:
            subdir          = file.split("/")[-1]
            subFromPath     = os.path.join( file,   "*" )
            subToPath       = os.path.join( toPath, subdir+"/" )
#            logger.info( "Creating subdirectory %s"%subToPath )
            copyDPMFiles( subFromPath, subToPath, toLocal=toLocal, fromLocal=fromLocal )


if __name__ == "__main__":

    # User specific
    from Analysis.Tools.user import dpm_directory

    def get_parser():
        ''' Argument parser for post-processing module.
        '''
        import argparse
        argParser = argparse.ArgumentParser(description = "Argument parser for nanoPostProcessing")
        argParser.add_argument('--ls',            action='store',      type=str,          help="list dpm path content")
        argParser.add_argument('--cp',            action='store',      type=str, nargs=2, help="copy dpm files")
        argParser.add_argument('--fromLocal',     action='store_true',                    help="copy from local directory")
        argParser.add_argument('--toLocal',       action='store_true',                    help="copy to local directory")
        argParser.add_argument('--mkdir',         action='store',      type=str,          help="make dpm dir")
        argParser.add_argument('--rm',            action='store',      type=str,          help="dpm path to remove")
        argParser.add_argument('--checkRootFile', action='store_true',                    help="check root file if corrupt")
        return argParser

    args = get_parser().parse_args()

    # Logging
    if __name__=="__main__":
        import Analysis.Tools.logger as logger
        logger = logger.get_logger("INFO", logFile = None )
        import RootTools.core.logger as logger_rt
        logger_rt = logger_rt.get_logger("INFO", logFile = None )

    else:
        import logging
        logger = logging.getLogger(__name__)


    if args.ls:
        if not args.ls.startswith("/dpm/"):
            args.ls = os.path.join( dpm_directory, args.ls )
        if args.ls.endswith("/"): args.ls += "*"
        allDirs = convertToPathList( args.ls )
        for dir in allDirs:        
            print dir

    elif args.cp:
        fromPath, toPath = args.cp

        if not fromPath.startswith("/dpm/") and not args.fromLocal:
            fromPath = os.path.join( dpm_directory, fromPath )
        if not toPath.startswith("/dpm/") and not args.toLocal:
            toPath = os.path.join( dpm_directory, toPath )

        copyDPMFiles( fromPath, toPath, toLocal=args.toLocal, fromLocal=args.fromLocal, rootFileCheck=args.checkRootFile )

    elif args.mkdir:
        if not args.mkdir.startswith("/dpm/"):
            args.mkdir = os.path.join( dpm_directory, args.mkdir )

        makeDPMDir( args.mkdir )

    elif args.rm:
        if not args.rm.startswith("/dpm/"):
            args.rm = os.path.join( dpm_directory, args.rm )

        removeDPMFiles( args.rm )

