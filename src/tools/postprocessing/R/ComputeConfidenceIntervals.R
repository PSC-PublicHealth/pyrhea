#!/usr/bin/env Rscript

library(pacman)
p_load(argparse,proto,findpython,getopt,jsonlite,tidyverse,data.table,boot,foreach,doMC,lubridate)

parser = ArgumentParser(description='Compute Confidence Intervals')

parser$add_argument('--input_directory', dest='input.directory', action='store', required=T)
parser$add_argument('--output_directory', dest='output.directory', action='store', required=T)
parser$add_argument('--scenario_names', dest='scenario.names', action='store', required=T)
parser$add_argument('--column_names', dest='column.names', action='store', required=F,
                    default=paste(c('Prev within 13','Prev outside 13','Prev within Cook','Prev outside Cook',
                    'Prev target','Prev nonTarget','Prev regionWide','Inc within 13','Inc outside 13',
                    'Inc within Cook','Inc outside Cook','Inc target','Inc nonTarget','Inc regionWide'), collapse=','))
parser$add_argument('--tier_value_names', dest='tier.value.names', action='store', required=F,
                    default=paste(c(
                      #  'Colonized Patient Days',
                      #  'Patient Bed Days',
                      #  'CRE Baths Given Out',
                      #  'CRE Swabs Used',
                      #  'Number of Patients Put on CP',
                      #  'CRE CP Days due to passive surveillance',
                      #  'CRE CP Days due to active surveillance',
                      #  'CRE CP Days due to xdro registry',
                      #  'CP Days for other reasons',
                      #  'XDRO Admissions',
                      'Newly Colonized',
                      'CRE Colonized Patients Admissions',
                      'Patient Admissions',
                      'Contact Precaution Days'
                    ), collapse=','))
parser$add_argument('--prevalence_and_incidence_pattern', dest='prevalence.and.incidence.pattern', action='store', required=F,
                    default='raw_prevalence_and_incidence_pre_day_13mile')
parser$add_argument('--raw_by_tier_pattern', dest='raw.by.tier.pattern', action='store', required=F,
                    default='raw_by_tier')

args = parser$parse_args()

input.directory = args$input.directory
output.directory = args$output.directory
dir.create(output.directory)

scenario.names = unlist(strsplit(args$scenario.names, ','))
column.names = unlist(strsplit(args$column.names, ','))
tier.value.cols = unlist(strsplit(args$tier.value.names, ','))

raw.revalence.and.incidence.pattern = args$prevalence.and.incidence.pattern
raw.by.tier.pattern = args$raw.by.tier.pattern



#scenario.names = c(
#  "Baseline_CRE_Prevalence_13_Mile_cumulative_counts", "Baseline_CRE_Prevalence_No_Constraint_cumulative_counts",
#  "Baseline_LTACH_vSNF_13_Mile_cumulative_counts", "Baseline_LTACH_vSNF_No_Constraints_cumulative_counts",
#  "Baseline_Largest_13-Mile_cumulative_counts", "Baseline_Largest_No_Constraints_cumulative_counts",
#  "Baseline_Min_LTAC_CRE_Prev_13_Mile_cumulative_counts", "Baseline_Min_SNA_CRE_13_Mile_cumulative_counts",
#  "Baseline_Random_13_Mile_cumulative_counts", "Baseline_Random_No_Constraints_cumulative_counts",
#  "Baseline_SNA_CRE_Prev_13_Mile_cumulative_counts", "Baseline_SNA_Power_13-Mile_cumulative_counts",
#  "Baseline_SNA_Power_No_Constraints_cumulative_counts", "CRE_Prevalence_13_Mile_cumulative_counts",
#  "CRE_Prevalence_No_Constraint_cumulative_counts", "LTACH_vSNF_13_Mile_cumulative_counts",
#  "LTACH_vSNF_No_Constraints_cumulative_counts", "Largest_13-Mile_cumulative_counts",
#  "Largest_No_Constraints_cumulative_counts", "Min_LTAC_CRE_Prev_13_Mile_cumulative_counts",
#  "Min_SNA_CRE_13_Mile_cumulative_counts", "Random_13_Mile_cumulative_counts", "Random_No_Constraints_cumulative_counts",
#  "SNA_CRE_Prev_13_Mile_cumulative_counts", "SNA_Power_13-Mile_cumulative_counts", "SNA_Power_No_Constraints_cumulative_counts")


read.data = function(directory, file.type) {
  foreach(scenario=scenario.names, .combine=rbind) %do% {
    d = fread(paste0(directory,'/',scenario,'_',file.type,'.csv'))
    d$scenario = scenario
    d
  }
}

value.cols = list(Names=column.names)
value.cols$outside = regexpr('outside', value.cols$Names) > 0
value.cols$Prev = regexpr('Prev', value.cols$Names) > 0

starting.date = ymd('2013-01-01')

d.grouped = read.data(input.directory, raw.revalence.and.incidence.pattern)
min.Day = min(d.grouped$Day)
d.grouped[,Day.Interval.30:=1+as.integer((Day-min.Day)/30)]
d.grouped[,Date:=starting.date+Day]

registerDoMC(8)

melt(d.grouped[, c('scenario','run','Day.Interval.30', value.cols$Names), with=F],
     id.vars=c('scenario','run','Day.Interval.30'), measure.vars = value.cols$Names) -> d.grouped.long

grouped.data.stats = function(d) {
  b = boot(d$value, R=nrow(d), statistic=function(x,i) median(x[i]))
  b.ci = boot.ci(b, conf=0.95, type='basic')
  
  d.by.run = d[, .(mean.value=mean(value)), by='run']
  mean.by.run = as.list(summary(d.by.run$mean.value))
  mean.by.run$IQR = IQR(d.by.run$mean.value)
  names(mean.by.run) = paste0('mean.by.run.', names(mean.by.run))
   
  ungrouped.daily = as.list(summary(d$value))
  ungrouped.daily$IQR = IQR(d$value)
  names(ungrouped.daily) = paste0('ungrouped.daily.', names(ungrouped.daily))
  
  c(list(
      ungrouped.daily.basic.bootstrap.median = b.ci$t0,
      ungrouped.daily.basic.bootstrap.95.ci.lower = b.ci$basic[4],
      ungrouped.daily.basic.bootstrap.95.ci.upper = b.ci$basic[5],
      ungrouped.basic.bootstrap.95.ci.width = b.ci$basic[5] - b.ci$basic[4]),
    ungrouped.daily,
    mean.by.run
  )
}

foreach(scenario.name=scenario.names, .combine=rbind) %dopar% {
  d.grouped.long[scenario==scenario.name,
                 grouped.data.stats(.SD),
                 by=.(scenario, variable, Day.Interval.30),
                 .SDcols=c('run','value')]
} -> d.final


d.tier = read.data(input.directory, raw.by.tier.pattern)

melt(d.tier[, c('scenario','run','Tier of Care', tier.value.cols), with=F],
     id.vars=c('scenario','run','Tier of Care'), measure.vars = tier.value.cols) -> d.tier.long

foreach(scenario.name=scenario.names, .combine=rbind) %dopar% {
  d.tier.long[scenario==scenario.name,
                 grouped.data.stats(.SD),
                 by=c('scenario', 'variable', 'Tier of Care'),
                 .SDcols=c('run','value')]
} -> d.tier.final


d.extra = d.grouped[,
                    .(Inc.Within=sum(`Inc within Cook`), Inc.Region.Wide=sum(`Inc within Cook`+`Inc outside Cook`)),
                    by=c('scenario','run')]

d.extra.long = melt(d.extra, id.vars=c('scenario','run'))

foreach(scenario.name=scenario.names, .combine=rbind) %dopar% {
  d.extra.long[scenario==scenario.name,
              grouped.data.stats(.SD),
              by=c('scenario','variable'),
              .SDcols=c('run','value')]
} -> d.extra.final


fwrite(d.final, paste(output.directory, 'prev.and.inc.csv', sep=','))
fwrite(d.tier.final, paste(output.directory, 'by.tier.csv', sep='/'))
fwrite(d.extra.final, paste(output.directory, 'incidence.within.cook.and.regionwide.csv', sep='/'))
    



