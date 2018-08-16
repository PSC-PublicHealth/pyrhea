#! /usr/bin/bash -x

topdir='/pylon5/pscstaff/welling/pyrhea/cre_bundle_201808'
outbasedir='/pylon5/pscstaff/jleonard/pyrhea'
pyrheadir='/pylon5/pscstaff/welling/git/pyrhea'
runyaml='ChicagoLand_xdrobundle_5_3.yaml'
DEFAULT='scenario_names.txt'
minNotesSz=900000000

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
	| sed "s*%runyaml%*$runyaml*g"
}

scenariolist=${1-$DEFAULT}
for scenario in `cat $scenariolist`
do
    bundledir=$outbasedir/output/cre_bundle
    targetdir=$bundledir/baseline
    #targetdir=$bundledir/initial_scenarios/$scenario
    nNotes=`$topdir/find_notes.py $targetdir $minNotesSz | wc -l`
    maxNote=$(( $nNotes-1 ))

    mkdir -p baseline_$scenario
    pushd baseline_$scenario > /dev/null
    echo 'building ' $PWD
    $topdir/gen_facility_yaml.py $bundledir/$runyaml $scenario \
	> xdro_facs.yaml
    customize $topdir/gen_counts_parallel.proto > gen_counts_parallel.sl
    customize $topdir/gen_csvfiles.proto > gen_csvfiles_parallel.sl
    customize $topdir/gen_costs_parallel.proto > gen_costs_parallel.sl
    echo *
    popd > /dev/null
done

