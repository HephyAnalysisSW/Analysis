''' Implementation of a directory based results DB for CMS analyses
    Supports merging and does not destroy afs volumes.

'''

# Standard imports
import os
import pickle
import uuid

# Logger
import logging
logger = logging.getLogger(__name__)

# Try to read a result from a single file. 
def read_from_file( f, key = None, forgiving = True):
    if os.path.exists( f ):
        try:
            with open(f) as _f:
                res = pickle.load(_f)
            if key is not None:
                if res.has_key( key ):
                    return res[ key ]
            else:
                return res
        except IOError:# Nothing found
            logger.warning( "Warning! Ignoring IOError when reading %s", f)
            pass 
        except ValueError:# It can be that we're loading from a tmp file that's currently being written. This gives 'ValueError: insecure string pickle'
            logger.warning( "Warning! Ignoring ValueError when reading %s", f)
            pass
        except EOFError:#same
            logger.warning( "Warning! Ignoring EOFError when reading %s", f)
            pass 
        except Exception as e: #something else wrong?
            logger.error( "Error reading file %s", f )
            if not forgiving:
                raise e
    return None

class MergingDirDB:
    def __init__( self, directory, init_on_start = True):
        '''
        Will create the directory if it doesn't exist 
        '''
        # work directory
        self.directory = directory
        # dictinary where the instance stores its unique 

        # create directory
        if not os.path.isdir( self.directory ):
            # there can be a race condition, if many jobs starts at the same time and also create this 
            # directory
            try:
                os.makedirs( self.directory ) 
            except OSError as e:
                if not os.path.isdir( self.directory):
                    raise e

        # new data
        self.data_dict = {}
        # read all files when starting?
        if init_on_start:
            self.data_on_disk_dict = self.data_from_all_files()
        else:
            self.data_on_disk_dict = {}

        # make a unique pkl file to write for this process
        self.unique_tmp_file = 'tmp_'+str(uuid.uuid4()) 

## Let's try WITHOUT done/tmp for the moment
#    def flush( self ):
#        # Rename the current tmp_file to 'done_...' and make a new tmp file
#        tmp_file  = os.path.join( self.directory, self.unique_tmp_file )
#        
#        done_file = os.path.join( self.directory, 'done_' + self.unique_tmp_file.lstrip('tmp_') )
#
#        if os.path.exists( tmp_file ):
#            os.rename( tmp_file, done_file )
#        
#        self.unique_tmp_file = 'tmp_'+str(uuid.uuid4()) 

    def data_from_all_files( self ):
        '''Read from all files in increasing order of unix time.'''
        files = self.tmp_files()
        if os.path.exists( self.merged_file() ):
            files.append( self.merged_file() )
        files = [ (f, os.path.getmtime(f)) for f in files ]
        files.sort( key = lambda r:r[1] )
        data = {} 
        for f, _ in files:
            _data = read_from_file( f )
            if _data:
                data.update(_data)
                
        return data

    def add(self, key, data, overwrite=False):

        if not overwrite:
            if self.read_from_all_files(key) is not None:
                logger.warning( "Already found key %r . Do not store data.", key )
                return data

        # Add data to private dictinary and store the file
        self.data_dict[key] = data 
        try:
            with open(os.path.join( self.directory, self.unique_tmp_file), 'w') as _f:
                pickle.dump( self.data_dict, _f )
        except Exception as e:
            logger.error( "Something wrong with file %s",  os.path.join( self.directory, self.unique_tmp_file) )
            raise e
        logger.debug( "Added key %r to file %s", key, os.path.join( self.directory, self.unique_tmp_file) )
        return data

    def get(self, key):
        ''' Get all entries in the database matching the provided key.
        '''

        if self.data_on_disk_dict.has_key( key ): return self.data_on_disk_dict[ key ]
        if self.data_dict.has_key( key ): return self.data_dict[ key ]
        # if we don't alread have the key, load it from all files and remember it in case you're asked again:
        return self.read_from_all_files(key)

    def contains(self, key):
        ''' Get all entries in the database matching the provided key.
        '''
        return self.get(key) is not None

    # Here we collect all files from other processes. 
    def tmp_files( self ):
        return [ os.path.join( self.directory, f ) for f in os.listdir(self.directory) if f.startswith( 'tmp_') ]

    # Here we find  all files that are already merged.
    # Merging should happen while NO jobs are running (because they are reading the merge file) 
    def merged_file( self ):
        merged_file = os.path.join( self.directory, 'merged' )
        return merged_file

    def read_from_all_files( self, key ):
        '''Read from all files that could possibly contain
           the result and return the newest according to unix modification time'''
        results = [] 
        for f in self.tmp_files() + [self.merged_file()]:
            result =  read_from_file( f, key )
            if result is not None:
                results.append( [ result, os.path.getmtime( f ) ] )
        if len(results)==0: return None
        results.sort( key = lambda r:r[1] )
        if len(results)>1:
            logger.warning( "Found %i results with different timestamp for key %r . Return newest.", len(results), key )
        return results[-1][0]
        
    def merge( self, clear = False):
        if len(self.tmp_files())==0:
            logger.info( "No tmp files, nothing to do.")
            return
        if os.path.exists( self.merged_file() ):
            f = self.merged_file()
            result = None
            if os.path.exists(f):
                try:
                    with open(f) as _f:
                        result = pickle.load(_f)
                except IOError:
                    pass
                except Exception as e:
                    logger.error( "Error reading file %s", f )
                    raise e
        else:
            result = {}
        logger.info( 'Found %i keys in merged file.', len(result.keys()) )
        # result will be 'None' if loading from an existing merge file failed.
        # That will result in an error below which is what we want, because in that case we don't want to write anything
        results = []
        for f in self.tmp_files():
            results.append( (f, os.path.getmtime( f )) )
        results.sort( key = lambda r:r[1] )
        for _result, _ in results:
            try:
                with open(_result) as _f:
                    result.update( pickle.load(_f) )
            except Exception as e:
                logger.error( "Something wrong with file %s", _result)
                raise e
        try:
            with open(self.merged_file(), 'w') as _f:
                pickle.dump(result, _f)
        except Exception as e:
            logger.error( "Something wrong with file %s", self.merged_file() )
            raise e
        logger.info( 'Wrote %i keys to merged file.', len(result.keys()) )
        
        if clear and os.path.exists( self.merged_file() ):
            try:
                with open(self.merged_file()) as _f:
                    pickle.load(_f)
            except:
                logger.error( "Could not load merged pickle file %s. Will not delete tmp files.",  self.merged_file() )
                return
            logger.info( "Merged pkl file %s seems OK, will delete tmp files.", self.merged_file() )
            for f in self.tmp_files():
                os.remove( f )
        return True

if __name__ == "__main__":
    import Analysis.Tools.logger as logger
    logger    = logger.get_logger( "DEBUG", logFile = None)

    import ROOT

    dirDB = MergingDirDB("./test2")

    #dirDB.add('x',1)
    #dirDB.add('y',2, overwrite=True)
    #dirDB.merge()
