''' Runs delphes on an edm root file
See https://twiki.cern.ch/twiki/bin/viewauth/CMS/DelphesUPG
'''

# Standard imports
import os
import subprocess
import ROOT

# Logger
import logging
logger = logging.getLogger(__name__)

# cmd line syntax:
#./DelphesCMSFWLite cards/delphes_card_CMS.tcl out.root /afs/hephy.at/data/rschoefbeck02/TTXPheno/small/ttZ0j_rwgt_patch_625_highStat_slc6_amd64_gcc630_CMSSW_9_3_0_tarball_small.root

class DelphesProducer:

    def __init__(self, card = 'delphes_card_CMS'):
        self.delphes_dir = os.path.expandvars( '$CMSSW_BASE/../delphes' )
        self.card = 'cards/' + card + '.tcl'

    def produce( self, infiles, outfile, executable = 'fwlite'):


        if executable == 'hepmc':
            executable_ = 'DelphesHepMC'
        elif executable == 'fwlite':
            executable_ = 'DelphesCMSFWLite'
        
        logger.info( "Running %s on %i files using %s", executable_, len(infiles), self.card )

        # Clean output
        tmp_files = []
        if os.path.exists( outfile ):
            logger.warning( "Found output file %s. Deleting", outfile )
            os.remove( outfile )

        # Create output directory if necessary
        if not os.path.exists( os.path.dirname(outfile) ):
            try:
                os.makedirs(os.path.dirname(outfile))
            except OSError:
                pass

        # process files individually
        for i_infile, infile in enumerate(infiles):
            tmp_file = outfile.rstrip('.root')+'_tmp_%i.root'%i_infile
            if os.path.exists( tmp_file ):
                logger.warning( "Found output file %s. Deleting", tmp_file )
                os.remove( tmp_file )
            tmp_files.append( tmp_file )
            logger.info( "Calling Delphes: %s", " ".join( [ os.path.join( self.delphes_dir, executable_), os.path.join( self.delphes_dir, self.card), tmp_file, infile ] ) )
            subprocess.check_call( [ os.path.join( self.delphes_dir, executable_), os.path.join( self.delphes_dir, self.card ), tmp_file, infile ] )

        # Hadd
        subprocess.check_call( ['hadd', outfile] + tmp_files ) 

        # Clean up 
        for tmp_file in tmp_files:
            os.remove(tmp_file)

