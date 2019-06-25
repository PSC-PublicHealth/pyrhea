# Notes on Setting Up and Running RHEA on bridges.psc.edu #

## To set up passwordless ssh transfers from bridges to github ##
```
cd ~/.ssh
ssh-keygen -o
```
Then input the key into github online

## To set up a working directory on /pylon5 ##

The home directory is on a small NFS filesystem, so it's useful to make a stable workspace on the /pylon5 filesystem.  If your username happens to be lwinch and your group happens to be ca5phmp, do:
```
ln -s /pylon5/ca5phmp/lwinch ~/WORK # create virtual link if you are working on supercomputer
cd ~/WORK
mkdir rhea
cd rhea
git clone --recurse-submodules https://github.com/PSC-PublicHealth/pyrhea.git
git clone https://github.com/PSC-PublicHealth/quilt.git
git clone https://github.com/PSC-PublicHealth/phacsl.git
```

## Work through creating a python environment suitable for running RHEA ##
1. Create a pyrhea python env:
```
module load anaconda5
conda create -n pyrheaEnv python=2.7 scipy numpy pip yaml six matplotlib
```
2. Create a file that will be sourced to bring the environment into effect
```
vim setup_python_path_olympus.sh
```
- make sure to include 'module load anaconda5' and 'source activate pyrheaEnv'
- make sure the python path has rhea sim, rhea tools, phacsl, and quilt
- for example, a working setup file might be:
```
#! /bin/bash
module load gcc/7.3.0
module load mpi/gcc_openmpi
module load anaconda5
source activate pyrheaEnv
ROOTDIR=/home/lwinch/WORK/rhea
export PYTHONPATH=${ROOTDIR}/pyrhea/src/sim:${ROOTDIR}/pyrhea/src/tools:${ROOTDIR}/phacsl/phacsl-utils/src:${ROOTDIR}/quilt/src
```
Once that file has been created, you can set up the environment with:
```
source setup_python_path_olympus.sh 
```
You dont need to this normally because the slurm script will try to do it automatically

3. try to start running the following things to see if you are missing any libraries, and install them to the pyrhea env
```
cd pyrhea/src/sim/
./pyrhea.py --help
./pyrhea.py week_run_ChicagoLand.yaml
```
### alleged completed list: ###
scipy numpy pip yaml six matplotlib mpi4py mpich2 pika ujson pyyaml greenlet jsonschema python-lmdb msgpack bcolz retrying

## make a folder to organize experiment scenarios  ##
```
cd /pylon5/ca5phmp/lwinch/rhea
mkdir runs
cd runs
mkdir baseline_test
cd baseline_test
cp /pylon5/ca5phmp/lwinch/rhea/pyrhea/src/tools/postprocessing/* .
rm baseline_OC_prep.py  # you won't need that one
mv gen_scenario_custom_yaml_Chicago.py gen_scenario_custom_yaml.py
cp /pylon5/ca5phmp/welling/* .   # for a few files not yet in the github repo
```

There are some edits needed:
```
vim runs_by_array.proto
```
- make sure 'SBATCH --constraint="intel"' and 'SBATCH --mem-per-cpu=40000' are disabled. 
- choose job array range carefully.
- set 'SBATCH --ntasks-per-node=4'
- add space to fix syntax error line 35 char 66

The main script which customizes other scripts must be modified to the environment to be used for the runs.
```
vim generate_subdirs.bash
```
- edit with correct paths for pythonsetup and pyrheadir
- edit runyaml with correct intervention yaml file

```
vim scenario_names.txt
```
- make sure correct scenarios are uncommented

```
vim ChicagoLand_baseline_run_5_3.yaml
vim ChicagoLand_xdrobundle_5_3.yaml
```
- add correct modeldir

```
vim scenario_mod.py
```
- change the boolean at the top of the file that controls caching to local disk to True

That should be all the needed changes to files.

# The lockserver #
There is a separate program, lockserver.py, which helps to keep separately running parts of RHEA synchronized.  It's not needed for most operations, but it is needed if RHEA is to use local disks on the worker nodes or for taumod runs (for prevalence tuning).  
To start it up:
```
cd /pylon5/ca5phmp/lwinch/rhea/runs/baseline_test  # the directory containing generate_subdirs.bash
nohup /pylon5/ca5phmp/lwinch/rhea/pyrhea/src/sim/lockserver.py --debug > $workdir/lockserver.log 2>&1 &
```
This is done on the *front end*, rather than on an internal node- that way it can stay up more or less forever, and it takes very few resources.  This will cause lockserver.info to be created in that directory.

## Generate the Run Directories and Start Runs ##
```
bash generate_subdirs.bash
```
This reads lines from scenario_names.txt and creates corresponding subdirectories with appropriately customized files.

```
cd ltac3 # or whichever subdirs were created
sbatch runs_by_array.sl # to queue a job
squeue -u lwinch # to see your current jobs
```

## Defining New Facility Groups for Scenarios ##
```
vim pyrhea/models/ChicagoLand/scenario_constants/scenario_groups.csv 
```
- Edit here if you want to add facilities to specific scenarios or define new ones.  The column names are the scenario names.

## More Misc Operations ##

To see memory usage and better tailor jobs:
```
sacct -j jobid
```

To get an interactive session on one of the compute nodes (with more resources than the front end nodes):
```
interact
```

To see available allocations for your projects:
```
projects
```
