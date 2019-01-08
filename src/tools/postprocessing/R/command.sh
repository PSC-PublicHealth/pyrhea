#! /usr/bin/bash -x
scenario="$1"
root=${HOME}/sshfshook/201812_Chicago_xdro_bundle_post
#root=${HOME}/sshfshook/201811_Chicago_baseline_work2_post
echo 'scenario' $scenario

Rscript ComputeConfidenceIntervals.R \
	--input_directory $root/$scenario \
	--output_directory $root/$scenario \
	--scenario_names ${scenario}_cumulative_counts \
	--column_names 'Prev within 13,Prev outside 13,Prev within Cook,Prev outside Cook,Prev target,Prev nonTarget,Prev regionWide,Inc within 13,Inc outside 13,Inc within Cook,Inc outside Cook,Inc target,Inc nonTarget,Inc regionWide'
