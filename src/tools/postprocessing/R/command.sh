#! /usr/bin/bash -x
scenario="$1"
#root=${HOME}/sshfshook/201812_Chicago_xdro_bundle_post
#root=${HOME}/sshfshook/201811_Chicago_baseline_work2_post_early
#root=${HOME}/sshfshook/chicago_baseline_201903_post
#root=${HOME}/sshfshook/201811_Chicago_baseline_work2_redo_post
#root=${HOME}/sshfshook/chicago_xdro_201903
#root=${HOME}/sshfshook/chicago_crebundle_201904
#root=${HOME}/sshfshook/chicago_crebundle_baseline_201904_post
#root=${HOME}/sshfshook/201904_Chicago_ltac_cre_bundle
#root=${HOME}/sshfshook/chicago_xdrobundle_201906
#root=/tmp/csvfiles/201907_CRE_Bundle
root=${HOME}/sshfshook/ltach_scenarios_201907
echo 'scenario' $scenario

Rscript ComputeConfidenceIntervals.R \
	--input_directory $root/$scenario \
	--output_directory $root/$scenario \
	--scenario_names ${scenario}_cumulative_counts \
	--column_names 'Prev within 13,Prev outside 13,Prev within Cook,Prev outside Cook,Prev target,Prev nonTarget,Prev regionWide,Inc within 13,Inc outside 13,Inc within Cook,Inc outside Cook,Inc target,Inc nonTarget,Inc regionWide'
