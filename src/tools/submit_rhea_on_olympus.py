import csv
import sys,os
from optparse import OptionParser
import subprocess

def createPBSForRun(fileName,nreals,inputfile,outprefix,rheadir,jobName):
    
    with open(fileName,"wb") as f:
        f.write("#!/bin/bash\n")
        f.write("#PBS -l walltime=96:00:00\n")
        f.write("#PBS -l nodes=1:ppn=1\n")
        f.write("#PBS -q intel256\n")
        f.write("#PBS -l pmem=30gb\n")
        f.write("#PBS -t 1-{0}\n".format(nreals))
        f.write("#PBS -N {0}_run\n".format(jobName))
        
        f.write("module load io\n")
        f.write("module load openblas\n")
        f.write("module load python\n")
        f.write("source $WORK2/pyrhea-local/bin/activate\n")
        
        f.write("cd $PBS_O_WORKDIR\n")
        #f.write("sleep 100\n")
        f.write("python {0}/src/sim/pyrhea.py -o {1}_$PBS_ARRAYID.pkl {2} > out.{1}.$PBS_ARRAYID".format(rheadir,
                                                                                                         outprefix,
                                                                                                         inputfile))

def createPBSForAnalysis(fileName,runPBSId,inputfile,outprefix,rheadir,jobName):
    with open(fileName,"wb") as f:
        f.write("#!/bin/bash\n")
        f.write("#PBS -l walltime=2:00:00\n")
        f.write("#PBS -l nodes=1:ppn=56\n")
        f.write("#PBS -l mem=30gb\n")
        f.write("#PBS -q intel1024\n")
        f.write("#PBS -W depend=afterokarray:{0}[]\n".format(runPBSId))
        f.write("#PBS -N {0}_anal\n".format(jobName))
        f.write("module load io\n")
        f.write("module load openblas\n")
        f.write("module load python\n")
        f.write("source $WORK2/pyrhea-local/bin/activate\n")
        f.write("export OMP_NUMTHREADS=1\n")
        
        f.write("cd $PBS_O_WORKDIR\n")
        f.write("python {0}/src/tools/gatherCounts.py -o {2} -g -n '{2}*.pkl' -m 56 {1} > out.counts\n".format(rheadir,
                                                                                                              inputfile,
                                                                                                              outprefix,
                                                                                                          ))
        f.write("#python {0}/src/tools/compute_costs.py -o {2} -g -n '{2}*.pkl' -r 10000 -m 56 {1} > out.costs\n".format(rheadir,
                                                                                                                        inputfile,
                                                                                                                        outprefix,
                                                                                                                    ))
        
def main():
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl [--out outname.yaml] run_descr.yaml
    """)
    parser.add_option('-o', '--outprefix', action='store', type='string',
                     help="will be an output prefix for all output files",default="output_rhea")
    parser.add_option('-i', '--inputyaml', action='store', type='string',
                      help='input yaml file for this run')
    parser.add_option('-r', '--rheadir', action='store', type='string',
                      help='RHEA, RHEA, where the heck is RHEA',default="/mnt/lustre0/temporary/work/common/users/stbrown/pyrhea/")
    parser.add_option('-n', '--nreals', action='store', type='int',
                      help='Number of realizations to run Rhea for', default=20)
    parser.add_option('-j','--jobname',action='store',type='string',default="job")
    #parser.add_option('-t', '--tmpDir')
    
    opts, args = parser.parse_args()
    
    runFileName = "run_rhea_{0}.pbs".format(opts.inputyaml)
    
    createPBSForRun(runFileName, opts.nreals,opts.inputyaml,opts.outprefix,opts.rheadir,opts.jobname)
    cmd = 'qsub {0}'.format(runFileName)
    pbsId = subprocess.check_output([cmd],shell=True)
    
    ### Now that we have the pbsId, we can make the analysis pbs with a dependency

    print "pbsId = {0}".format(pbsId)
    sys.stdout.flush()
    analysisFileName = "run_anal_{0}.pbs".format(opts.inputyaml)
    print 
    createPBSForAnalysis(analysisFileName, pbsId, opts.inputyaml, opts.outprefix, opts.rheadir,opts.jobname)

    cmd = 'qsub {0}'.format(analysisFileName)
    aPBSId = subprocess.check_output([cmd],shell=True)
if __name__ == "__main__":
    main()

