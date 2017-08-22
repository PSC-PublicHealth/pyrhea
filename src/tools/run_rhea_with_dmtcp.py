import csv
import sys,os,shutil
import random
from optparse import OptionParser
import subprocess

port_start = 50000
port_end = 55000

def createBurninPBSFile(fileName,condaEnv,nreals,queue,inputfile,outprefix,dayCheck,
                        rheadir,jobName,checkpointDir,mainDir,memGB="20",walltime="96:00:00"):
     with open(fileName,"wb") as f:
        f.write("#!/bin/bash\n")
        f.write("#PBS -l walltime={0}\n".format(walltime))
        f.write("#PBS -l nodes=1:ppn=1\n")
        f.write("#PBS -q {0}\n".format(queue))
        f.write("#PBS -l mem={0}gb\n".format(memGB))
        f.write("#PBS -t 1-{0}\n".format(nreals))
        f.write("#PBS -N {0}_burnin_run\n".format(jobName))

        f.write("module load io\n")
        f.write("module load anaconda/4.3.1\n")
        f.write("module load dmtcp/trunk\n")
        f.write("source activate {0}\n".format(condaEnv))

        f.write("cd $PBS_O_WORKDIR\n")
        f.write("mkdir burnin_$PBS_ARRAYID; cd burnin_$PBS_ARRAYID\n")

        f.write("ln -sf {0}/{1} .\n".format(mainDir,inputfile))
        if os.path.exists("{0}/log_cfg.yaml".format(mainDir)):
            f.write("ln -sf {0}/log_cfg.yaml .\n".format(mainDir))
        if os.path.exists("{0}/cache".format(mainDir)):
            f.write("ln -sf {0}/cache .\n".format(mainDir))
        
        if os.path.exists("{0}/constants".format(mainDir)):
            f.write("ln -sf {0}/constants .\n".format(mainDir))
        
        f.write("dmtcp_launch --new-coordinator --rm python {0}/src/sim/pyrhea.py ".format(rheadir) 
                + "-k {0} -o {1}.pkl {2} > out.{1}.$PBS_ARRAYID\n".format(dayCheck,outprefix,inputfile))
        #f.write("ls ckpt_* | xargs -I {{}} ln -sf `pwd`/{{}} {0}/ \n".format(checkpointDir))
        #f.write("ls dmtcp_restart_script_* |  xargs -I {{}} ln -sf `pwd`/{{}} {0}/dmtcp_restart_script_$PBS_ARRAYID.sh\n".format(checkpointDir))
        f.write("checkjob $PBS_JOBID[$PBS_ARRAYID]\n")

def createRunPBSFile(fileName,condaEnv,nreals,nbreals,queue,
                     inputfile,outprefix,
                     rheadir,jobName,checkpointDir,resultsDir,mainDir,workDir,
                     pbsIdToHold,memGB="20",walltime="96:00:00"):

     with open(fileName,"wb") as f:
             f.write("#!/bin/bash\n")
             f.write("#PBS -l walltime={0}\n".format(walltime))
             f.write("#PBS -l nodes=1:ppn=1\n")
             f.write("#PBS -q {0}\n".format(queue))
             f.write("#PBS -l mem={0}gb\n".format(memGB))
             f.write("#PBS -t 1-{0}\n".format(nreals))
             f.write("#PBS -N {0}_run\n".format(jobName))
             f.write("#PBS -W depend=afterokarray:{0}\n".format(pbsIdToHold))

             f.write("module load io\n")
             f.write("module load anaconda/4.3.1\n")
             f.write("module load dmtcp/trunk\n")
             f.write("source activate {0}\n".format(condaEnv))

             f.write("cd $PBS_O_WORKDIR\n")
             f.write("mkdir run_$PBS_ARRAYID; cd run_$PBS_ARRAYID\n")
             f.write("mkdir tmp\n")
             f.write("mkdir tmp2\n")
             f.write("x=$(( ( RANDOM % {0}) + 1 ))\n".format(int(nbreals)))
	     f.write("echo Using $x\n")
	     f.write("p=$(( 50000 + $PBS_ARRAYID ))\n")
             f.write('pwd\n')
             f.write("dmtcp_coordinator --daemon --exit-on-last --tmpdir ./tmp2 --port $p \n")
             f.write("{0}/{1}/checkpoints/burnin_$x/dmtcp_restart_script.sh -p $p -h `hostname` --tmpdir ./tmp > out.log\n".format(mainDir,workDir))
             f.write("ln -sf `pwd`/{0}.pkl {2}/{1}/{0}_$PBS_ARRAYID.pkl\n".format(outprefix,resultsDir,mainDir))
             f.write("checkjob $PBS_JOBID[$PBS_ARRAYID]\n")

def createCountsPBSFile(fileName,condaEnv,queue,inputfile,outprefix,
                        rheadir, jobName, mainDir,
                        pbsIdToHold,memGB="200gb",walltime="96:00:00"):
     with open(fileName,"wb") as f:
             f.write("#!/bin/bash\n")
             f.write("#PBS -l walltime={0}\n".format(walltime))
             f.write("#PBS -l nodes=1:ppn=56\n")
             f.write("#PBS -q {0}\n".format(queue))
             f.write("#PBS -l mem={0}gb\n".format(memGB))
             f.write("#PBS -N {0}_run\n".format(jobName))
             f.write("#PBS -W depend=afterokarray:{0}\n".format(pbsIdToHold))

             f.write("module load io\n")
             f.write("module load anaconda/4.3.1\n")
             f.write("module load dmtcp/trunk\n")
             f.write("source activate {0}\n".format(condaEnv))

             f.write("cd {0}\n".format(mainDir))
             f.write("python {0}/src/tools/gatherCounts.py -o $PBS_O_WORDIR/{1} -g -n $PBS_O_WORKDIR/'*.pkl' -t -m 50 {2} > $PBS_O_WORKDIR/out.count".format(rheadir,outprefix,inputfile)
)
def main():  
 
    parser = OptionParser(usage="""
    %prog [-i inputyamlfile -n <number of run realizations> -b <number of burnin realizations> -k <day in which to checkpoint>
    """)
    parser.add_option('-o', '--outprefix', action='store', type='string',
                     help="will be an output prefix for all output files",default="output_rhea")
    parser.add_option('-i', '--input', action='store', type='string',
                      help='input yaml file for this run')
    parser.add_option('-r', '--rheadir', action='store', type='string',
                      help='RHEA, RHEA, where the heck is RHEA',default=None)
    parser.add_option('-n', '--nrunreals', action='store', type='int',
                      help='Number of realizations to run Rhea for', default=20)
    parser.add_option('-b', '--nburninreals', action='store', type='int',
                      help='Number of burning realizations to create, ideally way lower than --nreals')
    parser.add_option('-j','--jobname',action='store',type='string',default="job",
                      help='The job name prefix to place on all PBS submissionts')
    parser.add_option('-q','--queue',action='store',type='string',default='batch',
                      help='Specify a particular queue to use for these runs in PBS')
    parser.add_option('-k','--checkpointday',action='store',type='int',default=100,
                      help='Day in which to create the checkpoints to restart the runs')
    parser.add_option('-s','--nocleanup',action='store_true',default=False)
    parser.add_option('-f','--force',action='store_true',default=False)
    parser.add_option('-c','--condaenv',action='store',type='string',default='$WORK1/pyrhea-ana')


    opts,args = parser.parse_args()

    burninRunFileName = "run_rhea_burnin_{0}.pbs".format(opts.input)
    runFileName = "run_rhea_{0}.pbs".format(opts.input)
    countsFileName = "run_rhea_counts_{0}.pbs".format(opts.input)
    
    randomMon = "{0}".format(random.random())

    mainDir = os.getcwd()
    #workingDir = "temp.working.{0}".format(randomMon)
    workingDir = "{0}_workfiles".format(opts.outprefix)
    checkpointDir = "{0}_checkpoints".format(opts.outprefix)
    finalOutsDir = "{0}_outputs".format(opts.outprefix)
    

    for d in [workingDir,checkpointDir,finalOutsDir]:
        print "checking {0}".format(d)
        if os.path.exists(d):
            print "{0} exists".format(d)
            if not opts.force:
                raise RuntimeError("The Directory {0} already exists, if you want to override, rerun with -f".format(d))
            else:
                shutil.rmtree(d)
    
    for d in [workingDir,checkpointDir,finalOutsDir]:
        os.makedirs(d)
        
    ### create and submit the burnin jobs
    
    os.chdir("{0}/{1}".format(mainDir,workingDir))
    os.makedirs("checkpoints")
    os.chdir("checkpoints")
    
    createBurninPBSFile(burninRunFileName,opts.condaenv,opts.nburninreals,
                        opts.queue,opts.input,opts.outprefix,opts.checkpointday,
                        opts.rheadir,opts.jobname,
                        "{0}/{1}".format(mainDir,checkpointDir),mainDir)

    cmd = 'qsub {0}'.format(burninRunFileName)
    bPBSId = subprocess.check_output([cmd],shell=True)
    
    os.chdir("{0}/{1}".format(mainDir,workingDir))
    os.makedirs("runs")
    os.chdir("runs")

    createRunPBSFile(runFileName,opts.condaenv,opts.nrunreals,
                     opts.nburninreals,opts.queue,
                     opts.input,opts.outprefix,
                     opts.rheadir,opts.jobname,
                     checkpointDir,finalOutsDir,mainDir,workingDir,
                     bPBSId)

    cmd = 'qsub {0}'.format(runFileName)
    rPBSId = subprocess.check_output([cmd],shell=True)

    os.chdir("{0}/{1}".format(mainDir,finalOutsDir))
    createCountsPBSFile(countsFileName,opts.condaenv,opts.queue,opts.input,
                        opts.outprefix,
                        opts.rheadir, opts.jobname, mainDir,
                        rPBSId)
    
    cmd = 'qsub {0}'.format(countsFileName)
    cPBSId = subprocess.check_output([cmd],shell=True)

    
    

if __name__ == "__main__":
    main()
    

    
