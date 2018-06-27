library(pacman)
p_load(mixtools, data.table, tidyverse, lme4, MASS, fitdistrplus, car, foreach, doMC, plotly)

nh.meta = fread('nh-meta.csv', header=F, skip=1, col.names=c('name','code','admissions','beds','census'))
nh.meta = nh.meta[code!='FAIR']
nh.los = fread('nh-los.csv', header=F, skip=1, col.names=c('code','los','discharge.year'))
nh.los = nh.los[code!='FAIR']

# NOTE: we should drop SCNH and CPNH
nh.los = nh.los[code!='CPNH']
nh.los = nh.los[code!='SCNH']
# NOTE: we should fit EXTW separately
nh.los.extw = nh.los[code=='EXTW']
nh.los = nh.los[code!='EXTW']


nh.los[,los.uncorrected:=nh.los$los]
nh.los$los = nh.los$los - runif(nrow(nh.los))

censor.threshold = 1000
                                          
nh.los[,censored:=nh.los$los < censor.threshold]

init.lambda = c(0.239557,   0.760443)
init.shape =  c(0.615004,   1.367461)
init.scale =  c(154.628632, 20.287383)

fit = weibullRMM_SEM(nh.los$los, nh.los$censored, maxit=100, lambda=init.lambda, shape=init.shape, scale=init.scale, k=2, verb=T)

cat('summary of fitted mixture:')
summary(fit)

cat('deciles of corrected, quantized, observed data:')
quantile(ceiling(nh.los$los), seq(0,1,0.1))

cat('deciles of corrected, quantized, fitted data:')
quantile(ceiling(rweibullmix(n=nrow(nh.los), lambda = fit$lambda, shape = fit$shape, scale=fit$scale)), seq(0,1,0.1))

ggplot() + 
  geom_density(aes(x=nh.los$los),
               fill='gray',alpha='0.5') +
  geom_density(aes(x=rweibullmix(n=nrow(nh.los), lambda = fit$lambda, shape = fit$shape, scale=fit$scale)),
               alpha=0.5, color='green', fill='green') +
  scale_x_log10() + ggtitle('Log-scale Density of Observed (gray) and fitted (green) NH LOS') +
  xlab('Log-scale NH LOS')

mean.component.1 = round(mean(rweibull(10000000,shape=fit$shape[1],scale=fit$scale[1])),1)
mean.component.2 = round(mean(rweibull(10000000,shape=fit$shape[2],scale=fit$scale[2])),1)

d_ = nh.los[,.(mean.los=mean(los), median.los=median(los), n.obs=length(los), mean.component.1, mean.component.2, k=(mean(los)-mean.component.2)/(mean.component.1-mean.component.2)), by='code']
d_[k<0,k:=0]
d_[k>1,k:=1]
fwrite(d_,'per-facility-k.csv')
sink('dist-params.csv')
summary(fit)
sink()

qqdata = function(a, b, quantile.width=0.1) {
  data.frame(
    q=seq(0,1,quantile.width),
    a=melt(quantile(a, seq(0,1,quantile.width)))$value,
    b=melt(quantile(b, seq(0,1,quantile.width)))$value
  )
}

compare.by.facility = function(facility.code) {
  params.fac = as.list(d_[code==facility.code])
  d.fac = nh.los[code==facility.code]
  ggplot() + 
    geom_histogram(aes(x=rweibullmix(n=nrow(d.fac), lambda = params.fac$k, shape = fit$shape, scale=fit$scale)),
                 alpha=0.25, color='green', fill='green', binwidth=10) +
    geom_histogram(aes(x=d.fac$los),
                   fill='gray', alpha=0.5, binwidth=10) +
    #scale_x_log10() +
    ggtitle('Log-scale Density of Observed (gray) and fitted (green) NH LOS') +
    xlab('Log-scale NH LOS')
}

qqdata.by.facility = function(facility.code, w) {
  params.fac = as.list(d_[code==facility.code])
  d.fac = nh.los[code==facility.code]
  d.fac.fit = rweibullmix(n=nrow(d.fac), lambda = params.fac$k, shape = fit$shape, scale=fit$scale)
  qqdata(d.fac$los, d.fac.fit, w)
}
  
  
  