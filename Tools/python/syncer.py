# Logger
import logging
logger = logging.getLogger(__name__)

# module singleton to keep track of files
file_sync_storage = []

# Wrap TCanvas class Print function
import ROOT
class myTCanvas( ROOT.TCanvas ):
    # recall the argument
    def Print( self, *args):
        logger.debug( "Appending file %s", args[0] )
        file_sync_storage.append( args[0] )
         
        super(myTCanvas, self).Print(*args)
ROOT.TCanvas = myTCanvas 

# pickle dump
import pickle

pickle._dump = pickle.dump
def syncer_pickle_dump( *args ):
    # first argument is file handle!
    if len(args)>1:
        file_sync_storage.append( args[1].name )
    else:
        logger.warning( "Pickle dump called with less than two arguments... shouldn't happen." )
    pickle._dump(*args)

pickle.dump = syncer_pickle_dump
     

# What happens on exit 
def sync_files():
        # No logger here, since already unloaded!
        # write outfile.sh because we can't scp in the container with kerberos authentication
    scp_cmd_filename = 'file_sync_storage.txt'
    with file( scp_cmd_filename, 'w' ) as outfile:
        for filename in file_sync_storage:
            if 'www' in filename:
                outfile.write('{filename}\n'.format(filename=filename))
            else:
                print "Will not sync %s" % filename
    print "Written %s" % scp_cmd_filename 

import atexit
atexit.register( sync_files )
