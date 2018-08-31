library(pacman)
p_load(mixtools, data.table, tidyverse, lme4, MASS, fitdistrplus, car, foreach, doMC, plotly)

# Load metafiles.  nh-meta.csv corresponds to
# models/ChicagoLand/OC_Hospital_and_Nursing_Home_Characteristics_for_RHEA_2.0_-_Adult_Only_-_09-15-2017_FINAL_NH_CHAR.csv
# and nh-los.csv to 
# models/ChicagoLand/OC_Nursing_Home_LOS_Line-Lists_for_RHEA_2.0_-_2011-2015_-_Adult_Only_-_09-29-2017_FINAL_NH_LOS_Line_List.csv
# with minor character clean-up

nh.meta = fread('nh-meta.csv', header=F, skip=1, col.names=c('name','code','admissions','beds','census'))
nh.los = fread('nh-los.csv', header=F, skip=1, col.names=c('code','los','discharge.year'))

# Exclude excluded sites
excludeList = list("CPNH", "FAIR", "SCNH", "SJNH")
for (eFac in excludeList) {
	nh.meta = nh.meta[code != eFac]
	nh.los = nh.los[code != eFac]
}

nh.los[,los.uncorrected:=nh.los$los]
nh.los$los = nh.los$los - runif(nrow(nh.los))

censor.threshold = 1000
                                          
nh.los[,censored:=nh.los$los < censor.threshold]

init.lambda = c(0.239557,   0.760443)
init.shape =  c(0.615004,   1.367461)
init.scale =  c(154.628632, 20.287383)

# Separate out EXTW, which has separate characteristics
extw.meta = nh.meta[code == 'EXTW']
nh.meta = nh.meta[code != 'EXTW']
extw.los = nh.los[code == 'EXTW']
nh.los = nh.los[code != 'EXTW']

cat('performing fitting for bulk of samples\n')
baseFit = weibullRMM_SEM(nh.los$los, nh.los$censored, maxit=100, lambda=init.lambda, shape=init.shape, scale=init.scale, k=2, verb=T)

cat('summary of fitted mixture for bulk of samples:')
summary(baseFit)

cat('deciles of corrected, quantized, observed data:')
quantile(ceiling(nh.los$los), seq(0,1,0.1))

cat('deciles of corrected, quantized, fitted data:')
quantile(ceiling(rweibullmix(n=nrow(nh.los), lambda = baseFit$lambda, shape = baseFit$shape, scale=baseFit$scale)), seq(0,1,0.1))

ggplot() + 
		geom_density(aes(x=nh.los$los),
				fill='gray',alpha='0.5') +
		geom_density(aes(x=rweibullmix(n=nrow(nh.los), lambda = baseFit$lambda, shape = baseFit$shape, scale=baseFit$scale)),
				alpha=0.5, color='green', fill='green') +
		scale_x_log10() + ggtitle('Log-scale Density of Observed (gray) and fitted (green) NH LOS Base Data') +
		xlab('Log-scale NH LOS')


cat('performing fitting for EXTW samples\n')
extwFit = weibullRMM_SEM(extw.los$los, extw.los$censored, maxit=100, lambda=init.lambda, shape=init.shape, scale=init.scale, k=2, verb=T)

cat('summary of fitted mixture for EXTW samples:')
summary(extwFit)

cat('deciles of corrected, quantized, observed EXTW data:')
quantile(ceiling(extw.los$los), seq(0,1,0.1))

cat('deciles of corrected, quantized, fitted EXTW data:')
quantile(ceiling(rweibullmix(n=nrow(extw.los), lambda = extwFit$lambda, shape = extwFit$shape, scale=extwFit$scale)), seq(0,1,0.1))

ggplot() + 
		geom_density(aes(x=extw.los$los),
				fill='gray',alpha='0.5') +
		geom_density(aes(x=rweibullmix(n=nrow(extw.los), lambda = extwFit$lambda, shape = extwFit$shape, scale=extwFit$scale)),
				alpha=0.5, color='green', fill='green') +
		scale_x_log10() + ggtitle('Log-scale Density of Observed (gray) and fitted (green) NH LOS EXTW Data') +
		xlab('Log-scale NH LOS')

cat('decomposing known means into base components\n')
mean.component.1 = round(mean(rweibull(10000000,shape=baseFit$shape[1],scale=baseFit$scale[1])),1)
mean.component.2 = round(mean(rweibull(10000000,shape=baseFit$shape[2],scale=baseFit$scale[2])),1)

d_ = nh.los[,.(mean.los=mean(los), median.los=median(los), n.obs=length(los), mean.component.1, mean.component.2, k=(mean(los)-mean.component.2)/(mean.component.1-mean.component.2)), by='code']
d_[k<0,k:=0]
d_[k>1,k:=1]
fwrite(d_,'per-facility-k.csv')

sink('dist-params.csv')
summary(baseFit)
sink()

sink('extw-dist-params.csv')
summary(extwFit)
sink()

qqdata = function(a, b, quantile.width=0.1) {
  data.frame(
    q=seq(0,1,quantile.width),
    a=melt(quantile(a, seq(0,1,quantile.width)))$value,
    b=melt(quantile(b, seq(0,1,quantile.width)))$value
  )
}
