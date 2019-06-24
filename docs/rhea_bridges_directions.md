# hermes.psc.edu/hangout
# ca5phmp group code for pylon5

cd ~/.ssh
ssh-keygen -o
# input key into github online

ln -s /pylon5/ca5phmp/lwinch ~/WORK # create virtual link if you are working on supercomputer
cd ~/WORK
mkdir rhea
cd rhea
git clone --recurse-submodules https://github.com/PSC-PublicHealth/pyrhea.git
git clone https://github.com/PSC-PublicHealth/quilt.git
git clone https://github.com/PSC-PublicHealth/phacsl.git

## BRIDGES ##
# create a pyrhea python env
module load anaconda5
conda create -n pyrheaEnv python=2.7 scipy numpy pip yaml six matplotlib
vim setup_python_path_olympus.sh
# make sure the python path has rhea sim, rhea tools, phacsl, and quilt
# make sure include 'module load anaconda5' and 'source activate' pyrhea env
source setup_python_path_olympus.sh # you dont need to this normally because the sl script will try to do it automatically

# try to start running the following things to see if you are missing any libraries, and install them to the pyrhea env
# alleged completed list: scipy numpy pip yaml six matplotlib mpi4py mpich2 pika ujson pyyaml greenlet jsonschema python-lmdb msgpack bcolz retrying
cd pyrhea/src/sim/
./pyrhea.py --help
./pyrhea.py week_run_ChicagoLand.yaml

# make folder to organize experiment scenarios
cd /pylon5/ca5phmp/lwinch/rhea
mkdir runs
cd runs
mkdir baseline_test
cd baseline_test
cp /pylon5/ca5phmp/lwinch/rhea/pyrhea/src/tools/postprocessing/* .
rm baseline_OC_prep.py  # you won't need that one
mv gen_scenario_custom_yaml_Chicago.py gen_scenario_custom_yaml.py
cp /pylon5/ca5phmp/welling/* .

vim runs_by_array.proto
# make sure 'SBATCH --constraint="intel"' and 'SBATCH --mem-per-cpu=40000' are disabled. 
# choose job array range carefully.
# set 'SBATCH --ntasks-per-node=4'
# add space to fix syntax error line 35 char 66

vim generate_subdirs.bash
# edit with correct paths for pythonsetup and pyrheadir
# edit runyaml with correct intervention yaml file

vim scenario_names.txt
# make sure correct scenarios are uncommented

vim ChicagoLand_baseline_run_5_3.yaml
vim ChicagoLand_xdrobundle_5_3.yaml
# add correct modeldir

# lockserver
1) Starting lockserver.py, like 'nice nohup python ./lockserver.py >& lockserver.log &' from your parent directory (corresponding to baseline_test rather than baseline_test/ltac3)
or nohup /pylon5/ca5phmp/lwinch/rhea/pyrhea/src/sim/lockserver.py --debug > $workdir/lockserver.log 2>&1 &
This is done on the *front end*, rather than on an internal node- that way it can stay up more or less forever, and it takes very few resources.  This will cause lockserver.info to be created in that directory.
2) Change the boolean in the top of scenario_mod.py that controls caching to local disk to True
That should be it.

bash generate_subdirs.bash

cd ltac3 # or whichever subdirs were created
sbatch runs_by_array.sl # to queue a job
squeue -u lwinch # to see your current jobs

vim pyrhea/models/ChicagoLand/scenario_constants/scenario_groups.csv 
# edit here if you want to add facilities to specific scenarios

# more misc operations
# sacct -j jobid to see mem usage and better taylor jobs
# interact to get more computing resources
# 'projects' to see allocation