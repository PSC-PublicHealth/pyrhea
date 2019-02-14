#! /usr/bin/bash -x

topdir='/pylon5/pscstaff/welling/pyrhea/oc_mrsa_201811b'
outbasedir=$topdir
pyrheadir='/pylon5/pscstaff/welling/git/pyrhea'
runyaml='OC_baseline_run_5_3.yaml'
prepscript='baseline_OC_prep.py'
DEFAULT='scenario_names.txt'
minNotesSz=100000000

customize () {
    cat $1 \
	| sed "s*%basedir%*$outbasedir*g" \
	| sed "s*%scenario%*$scenario*g" \
	| sed "s*%nfiles%*$nNotes*g" \
	| sed "s*%maxfile%*$maxNote*g" \
	| sed "s*%targetdir%*$targetdir*g" \
	| sed "s*%workdir%*${PWD}*g" \
	| sed "s*%bundledir%*$bundledir*g" \
	| sed "s*%pyrheadir%*$pyrheadir*g" \
	| sed "s*%topdir%*$topdir*g" \
	| sed "s*%minNotesSz%*$minNotesSz*g" \
	| sed "s*%runyaml%*$runyaml*g" \
	| sed "s*%prepscript%*$prepscript*g"
}

scenariolist=${1-$DEFAULT}
for scenario in `cat $scenariolist | grep -v '^#'`
do
    mkdir -p $scenario

    bundledir=$outbasedir
    #targetdir=$bundledir/baseline
    targetdir=$bundledir/$scenario
    nNotes=`$topdir/find_notes.py $targetdir $minNotesSz | wc -l`
    altNNotes=`$topdir/find_bcz.py $targetdir | wc -l`
    if [ $altNNotes > $nNotes ]; then
    	nNotes=$altNNotes
	fi
    if [ $nNotes -ne 0 ]; then
	maxNote=$(( $nNotes-1 ))
    else
	maxNote=50
    fi

    pushd $scenario > /dev/null
    echo 'building ' $PWD
    # $topdir/gen_facility_yaml.py $bundledir/$runyaml $scenario \
    #	> xdro_facs.yaml
    customize $topdir/runs_by_array.proto > runs_by_array.sl
    customize $topdir/gen_counts_parallel.proto > gen_counts_parallel.sl
    customize $topdir/gen_csvfiles.proto > gen_csvfiles_parallel.sl
    customize $topdir/gen_costs_parallel.proto > gen_costs_parallel.sl
    customize $topdir/run_diagnostics.proto > run_diagnostics.sl
    customize $topdir/cleanup.proto > run_diagnostics.sl
    customize $topdir/run_taumod.proto > run_taumod.sl
    customize $topdir/taumod_runner.proto > taumod_runner.bash
    customize $topdir/cleanup_taumod.proto > run_diagnostics.sl
    chmod +x taumod_runner.bash
    echo *
    popd > /dev/null
done

