#!/bin/sh
python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/example-ntuples-ttZ_dy --sample DY   --config_module Analysis.MVA.cfg_examples --config ttZ_dy_example --small
python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/example-ntuples-ttZ_dy --sample TTZ  --config_module Analysis.MVA.cfg_examples --config ttZ_dy_example --small
