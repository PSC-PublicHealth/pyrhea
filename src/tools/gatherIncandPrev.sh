#!/bin/bash

FILES=/home/jleonard/pyrhea/pyrhea/src/sim/runs/exp_scenlist_000009/*.pkl
for f in $FILES
do
  echo $f
  python ../tools/print_counts.py -n $f /home/jleonard/pyrhea/pyrhea/src/sim/runs/exp_scenlist_000009/year_allfac_ChicagoLand_crebundleXdro.yaml > $f.counts.out &  
done

wait
