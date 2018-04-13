#! /bin/bash -x

runid=$1
inrecs=$2

workdir=/tmp
tools=$HOME/git/pyRHEA_github/src/tools
#rundesc=$tools/../sim/twoyear_run_ChicagoLand.yaml
rundesc=$tools/../sim/twoyear_allfac_OC.yaml

myconvert () {
  low=$1
  high=$2
  $tools/arrival_time_parser.py -i $workdir/arrivals.txt -L $low -H $high \
    $rundesc
  mv arrival_time_arrays.npz arrays_indirect_${low}_${high}_${runid}.npz
}

grep DEBUG $inrecs | grep arrive > $workdir/arrivals.txt
#myconvert 0 100
myconvert 100 830
myconvert 100 465
myconvert 466 830


