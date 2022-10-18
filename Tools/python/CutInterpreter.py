''' Base class for cut interpreter 
'''

import logging
logger = logging.getLogger(__name__)

class CutInterpreter:
    ''' Translate var100to200-var2p etc.
    '''

    def __init__( self, continous_variables, discrete_variables, special_cuts):
        self.continous_variables = continous_variables
        self.discrete_variables  = discrete_variables
        self.special_cuts        = special_cuts

    def translate_cut_to_string( self, string ):

        # special cuts
        if string in list(self.special_cuts.keys()): return self.special_cuts[string]

        # continous Variables
        for var, tree_var in self.continous_variables:
            if string.startswith( var ):
                num_str = string[len( var ):].replace("to","To").split("To")
                upper = None
                lower = None
                if len(num_str)==2:
                    lower, upper = num_str
                elif len(num_str)==1:
                    lower = num_str[0]
                else:
                    raise ValueError( "Can't interpret string %s" % string )
                res_string = []
                if lower: res_string.append( tree_var+">="+lower )
                if upper: res_string.append( tree_var+"<"+upper )
                return "&&".join( res_string )

        # discrete Variables
        for var, tree_var in self.discrete_variables:
            logger.debug("Reading discrete cut %s as %s"%(var, tree_var))
            if string.startswith( var ):
                # So far no njet2To5
                if string[len( var ):].replace("to","To").count("To"):
                    raise NotImplementedError( "Can't interpret string with 'to' for discrete variable: %s. You just volunteered." % string )

                num_str = string[len( var ):]
                # logger.debug("Num string is %s"%(num_str))
                # var1p -> tree_var >= 1
                if num_str[-1] == 'p' and len(num_str)==2:
                    # logger.debug("Using cut string %s"%(tree_var+">="+num_str[0]))
                    return tree_var+">="+num_str[0]
                # var123->tree_var==1||tree_var==2||tree_var==3
                else:
                    vls = [ tree_var+"=="+c for c in num_str ]
                    if len(vls)==1:
                      # logger.debug("Using cut string %s"%vls[0])
                      return vls[0]
                    else:
                      # logger.debug("Using cut string %s"%'('+'||'.join(vls)+')')
                      return '('+'||'.join(vls)+')'
        raise ValueError( "Can't interpret string %s. All cuts %s" % (string,  ", ".join( [ c[0] for c in self.continous_variables + self.discrete_variables] +  list(self.special_cuts.keys()) ) ) )

    def cutString( self, cut, select = [""], ignore = [], photonEstimated=False):
        ''' Cutstring syntax: cut1-cut2-cut3
        '''
        if cut is None: return '1.'

        cuts = cut.split('-')
        # require selected
        cuts = [c for c in cuts if any( sel in c for sel in select )]
        # ignore
        cuts = [c for c in cuts if not any( ign in c for ign in ignore )]

        cutString = "&&".join( map( self.translate_cut_to_string, cuts ) )

        if photonEstimated:
          for var in ['met_pt','met_phi','metSig','dl_mt2ll','dl_mt2bb']:
            cutString = cutString.replace(var, var + '_photonEstimated')

        return cutString
    
    def cutList ( cut, select = [""], ignore = []):
        ''' Cutstring syntax: cut1-cut2-cut3
        '''
        if cut is None: return ['1.']

        cuts = cut.split('-')
        # require selected
        cuts = [c for c in cuts if any( sel in c for sel in select )]
        # ignore
        cuts = [c for c in cuts if not any( ign in c for ign in ignore )]
        return [ self.translate_cut_to_string(cut) for cut in cuts ] 
        #return  "&&".join( map( CutInterpreter.translate_cut_to_string, cuts ) )

