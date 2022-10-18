''' Implementation of a directory based results DB for CMS analyses


'''

# Standard imports
import os
import pickle


# Logger
import logging
logger = logging.getLogger(__name__)

class DirDB:
    def __init__( self, directory ):
        '''
        Will create the directory if it doesn't exist 
        '''
        self.directory = directory
        try: # errors can appear in parallel processing
            if not os.path.isdir( self.directory ):
                os.makedirs( self.directory )
        except:
            pass

    def __get_filename( self, key ):
        filename = str(hash(key))
        if len(filename)>4:
            return filename[:4] + '/' + filename[4:]
        else:
            return filename

    def get(self, key):
        ''' Get all entries in the database matching the provided key.
        '''

        result = None
        try:
            result = pickle.load( open(os.path.join( self.directory, self.__get_filename(key) )) ) 
        except IOError:
            # nothing found
            pass

        return result

    def contains(self, key):
        ''' We got that thing?
        '''
        return os.path.exists( os.path.join( self.directory, self.__get_filename(key) ) )

    def add(self, key, data, overwrite=False):

        filename = os.path.join( self.directory, self.__get_filename(key)) 
        try: # errors can appear in parallel processing
            if not os.path.isdir(os.path.dirname( filename )):
                os.makedirs( os.path.dirname( filename ) )
        except:
            pass

        if not overwrite:
            if os.path.exists( filename ):
                logger.warning( "Already found key '%r'. Do not store data.", key )
                return data
        with open( filename, 'wb' ) as _file:
            pickle.dump( data, _file )
        return data

if __name__ == "__main__":
    import Analysis.Tools.logger as logger
    logger    = logger.get_logger( "DEBUG", logFile = None)

    import ROOT

    dirDB = DirDB("./test")

    dirDB.add('y',1)
    dirDB.add(3,1)
    dirDB.add((2,3),ROOT.TH1F('x','x',100,0,1))
