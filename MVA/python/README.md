1. In your repository, create a directory <Repo/something/python/configs/>
   and put your MVA config there. An example is in Analysis/MVA/python/cfg_examples/ttZ_dy_example.py
   Also, copy the __init__.py  

2. Adapt the commands found in run_examples/make_ntuple.sh 
   to reflect the location of your config dir, the name of your config, output directory, and the names of the samples
   Run these commands to produce flat ntuples

3. Adapt the commands found in run_examples/multiclass.sh to perform the training 
