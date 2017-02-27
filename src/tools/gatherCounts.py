import csv
import sys,os
import glob
import numpy as np
import yaml
from multiprocessing import Process,Manager
from optparse import OptionParser

cwd = os.path.dirname(__file__)
sys.path.append(os.path.join(cwd, "../sim"))
import pyrheautils
import schemautils
import pathogenbase as pth
from facilitybase import CareTier
from notes_plotter import readFacFiles, checkInputFileSchema
from notes_plotter import SCHEMA_DIR, INPUT_SCHEMA
from time_series_plotter import mergeNotesFiles, getTimeSeriesList
import print_xdro_counts, print_counts

DEFAULT_OUT_FILE = 'counts_output.yaml'
def main():
    parser = OptionParser(usage="""
    %prog [--notes notes_file.pkl [--out outname.yaml] run_descr.yaml
    """)
    
    parser.add_option('-n', '--notes', action='append', type='string',
                     help="Notes filename - may be repeated")
    parser.add_option('-o', '--out', action='store', type='string',
                     help="Output file (defaults to %s)" % DEFAULT_OUT_FILE)
    parser.add_option('-g','--glob', action='store_true',
                      help=("Apply filename globbing for notes files."
                            "  (Remember to protect the filename string from the shell!)"))
    parser.add_option('-x','--xdroyamlfile',type='string', default=None,
                      help="specify the xdro scenario yaml file if there is one, if not, XDRO won't be costed")
    parser.add_option('-m','--nprocs',type='int',default=1,
                      help='number of cpus to run the costing model over')
    
    opts, args = parser.parse_args()
    
    inputDict = checkInputFileSchema(args[0], os.path.join(SCHEMA_DIR, INPUT_SCHEMA))
    if opts.out:
        outFileName = opts.out
    else:
        outFileName = DEFAULT_OUT_FILE
        
    

files = glob.glob('base1*.counts.out')
counts = {}
for file in files:
    with open(file,"rb") as f:
        csvReader = csv.reader(f)
	c = 0 
        for line in csvReader:
	    if c < 17:
		c += 1
		continue
            print line
            if len(line) > 1:
                abbrev = line[0]
                if abbrev not in counts.keys():
                   counts[abbrev] =[0,0]
                counts[abbrev][0] += float(line[2])
		counts[abbrev][1] += float(line[3])

averageCounts = {}
for abbrev,indprev in counts.items():
    print indprev
    indSum = indprev[0]/float(len(files)) 
    prevSum = indprev[1]/float(len(files)) 
    averageCounts[abbrev] = (indSum,prevSum)

with open('base2_average_counts.csv',"wb") as f:
    csvWriter = csv.writer(f)
    csvWriter.writerow(['Abbev','colonized (counts)','beddays (counts)'])
    for abbrev,values in averageCounts.items():
        csvWriter.writerow([abbrev,values[0],values[1]])
    
    


