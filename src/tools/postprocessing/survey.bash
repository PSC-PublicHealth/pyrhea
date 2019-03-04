#! /usr/bin/bash -x

#for scenario in `cat scenario_names.txt | grep -v '^#'`
for scenario in `cat scenario_names.txt | grep -v '^##' | sed 's/#//g'`
do
    echo $scenario
    echo '    mpz files:' `ls $scenario/*.mpz  2>/dev/null | wc -l`  
    csv_raw=`ls $scenario/*.csv 2>/dev/null | grep _raw_ | wc -l`
    csv_cost=`ls $scenario/*.csv 2>/dev/null | grep _costs_ | wc -l`
    csv_other=`ls $scenario/*.csv 2>/dev/null | grep -v _raw_ | grep -v _costs_ | wc -l`
    echo '    csv files: raw ' $csv_raw 'cost' $csv_cost 'other' $csv_other
    echo '    png files:' `ls $scenario/diagnostics/*.png 2>/dev/null | wc -l`
done
