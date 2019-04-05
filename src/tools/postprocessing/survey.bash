#! /usr/bin/bash -x

count_running() {
    jobid=$1
    echo $(squeue -u `whoami` | grep 'runs_' | grep "$jobid" | wc -l)
}

get_rhea_jobids() {
    if compgen -G "$scenario/slurm-*_*.out" > /dev/null; then
	echo $(for fname in $scenario/slurm-*_*.out; do s=`basename "$fname"`; s="${s##*-}"; s="${s%_*}"; echo "$s"; done | uniq)
    fi
}

for scenario in `cat scenario_names.txt | grep -v '^##' | sed 's/#//g'`
do
    echo '---------' $scenario
    . $scenario/run_info.bash
    lastday=$(( $totalrundays + 1 ))
    if compgen -G "$scenario/pyrhea_*_out.out"; then
	finished_runs=$(for fname in $scenario/pyrhea_*_out.out; do grep 'bump time' $fname | tail -1; done | grep $lastday | wc -l)
    else
	finished_runs=0
    fi
    num_running=0
    for jobid in $(get_rhea_jobids); do
	num_running=$(( $num_running + $(count_running $jobid) ))
    done
    echo '    RHEA runs: finished ' $finished_runs 'running' $num_running
    echo '    mpz files:' `ls $scenario/*.mpz  2>/dev/null | wc -l`  
    csv_raw=`ls $scenario/*.csv 2>/dev/null | grep _raw_ | wc -l`
    csv_cost=`ls $scenario/*.csv 2>/dev/null | grep _costs_ | wc -l`
    csv_other=`ls $scenario/*.csv 2>/dev/null | grep -v _raw_ | grep -v _costs_ | wc -l`
    echo '    csv files: raw ' $csv_raw 'cost' $csv_cost 'other' $csv_other
    echo '    png files:' `ls $scenario/diagnostics/*.png 2>/dev/null | wc -l`
done
