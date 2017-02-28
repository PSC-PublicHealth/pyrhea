import csv
import sys,os
from optparse import OptionParser

def createPBSForRun(fileName,nreals,inputfile,outprefix,rheadir):
    
    with open(fileName,"wb") as f:
        f.write("!/bin/bash\n")
        f.write("#PBS -l walltime=96:00:00\n")
        f.write("#PBS -l nodes=1:ppn=1\n")
        f.write("#PBS -q intel\n")
        f.write("#PBS -l pmem=30gb\n")
        f.write("#PBS -t 1-{0}\n".format(nreals))
        
        f.write("module load io\n")
        f.write("module load openblas\n")
        f.write("module load python\n")
        f.write("source $WORK2/pyrhea-local/bin/activate.csh\n")
        
        f.write("cd $PBS_O_WORKDIR\n")
        
        f.write("python {0}/src/sim/pyrhea.py -o {1}_$PBS_ARRAYID.pkl {2} > out.{1}.$PBS_ARRAYID".format(rheadir,
                                                                                                         outprefix,
                                                                                                         inputfile))

        
        
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
    #parser.add_option('-t', '--tmpDir')
    
    opts, args = parser.parse_args()
    
    runFileName = "run_rhea_{0}.pbs".format(opts.inputyaml)
    
    createPBSForRun(runFileName, opts.nreals,opts.inputyaml,opts.outprefix,opts.rheadir)
    
if __name__ == "__main__":
    main()

