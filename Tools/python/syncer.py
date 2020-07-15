'''
Wraps output functions and intercepts filenames.
All you have to do in your script is
import Analysis.Tools.syncer
and the filenames written via
    TCanvas::Print
    pickle.dump
go to 
    file_sync_storage.txt
which can then be synced.
'''

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
        # call original Print method 
        super(myTCanvas, self).Print(*args)
# what could possibly go wrong.
ROOT.TCanvas = myTCanvas 

# pickle dump
import pickle
# that's the old dump method
pickle._dump = pickle.dump
def syncer_pickle_dump( *args ):
    # second argument is file handle!
    if len(args)>1:
        file_sync_storage.append( args[1].name )
    else:
        logger.warning( "Pickle dump called with less than two arguments... shouldn't happen." )
    pickle._dump(*args)
# that's the new dump method
pickle.dump = syncer_pickle_dump

# What happens on exit 
def write_sync_files_txt():
    # No logger here, since it is already unloaded!
    # write outfile.sh because we can't scp in the container with kerberos authentication
    scp_cmd_filename = 'file_sync_storage.txt'
    if len(file_sync_storage)>0:
        with file( scp_cmd_filename, 'w' ) as outfile:
            for filename in file_sync_storage:
                if 'www/' in filename:
                    # for rsync cmd with relative path
                    outfile.write('{filename}\n'.format(filename=filename.replace('www/','www/./')))
                else:
                    print "Will not sync %s" % filename
        print "Analysis.Tools.syncer: Written %i files to %s for rsync." % (len(file_sync_storage), scp_cmd_filename)

import atexit
atexit.register( write_sync_files_txt )
