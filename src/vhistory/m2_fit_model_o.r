##############
## INITIATE ##
##############

source("f_betteRplots.r")
source("f_bngm.r")
source("f_model_o.r")
source("f_model_p.r")

library("Rcpp")
sourceCpp("f_demc.cpp")

#
###

##############
## INITIATE ##
##############

## Goal: load data, functions

## Time series
timeSeriesId = "all"

## Load data
MTS = read.table(paste("data/MTS_", timeSeriesId, ".csv", sep=""),sep=",",header=T)
head(MTS)

#
###

###########################
## FIT OBSERVATION MODEL ##
###########################

## Collectors
MTS_o = NULL
MTS_o_lo = NULL
MTS_o_hi = NULL
ddt.MTS_o = NULL
ddt.MTS_o_lo = NULL
ddt.MTS_o_hi = NULL
for (k in 1:length(unique(MTS[,1])))
{
  ## Time series index 
  idx = unique(MTS[,1])[k]
  
  ## Interpolate
  TS = MTS[MTS[,1] == idx,][,-1]
  list_o = fit_model_o(TS, sd1_o=.1, sd2_o=.1, K_o=30)
  
  ## Collect objects
  MTS_o = rbind(MTS_o, cbind(idx, list_o$Yhat_o))
  MTS_o_lo = rbind(MTS_o_lo, cbind(idx, list_o$Yhat_o_lo))
  MTS_o_hi = rbind(MTS_o_hi, cbind(idx, list_o$Yhat_o_hi))
  ddt.MTS_o = rbind(ddt.MTS_o, cbind(idx, list_o$ddt.Yhat_o))
  ddt.MTS_o_lo = rbind(ddt.MTS_o_lo, cbind(idx, list_o$ddt.Yhat_o_lo))
  ddt.MTS_o_hi = rbind(ddt.MTS_o_hi, cbind(idx, list_o$ddt.Yhat_o_hi))
}

# ## Save results
# system("mkdir out_o/")
# write.csv(x = MTS_o, file = "out_o/MTS_o.csv")
# write.csv(x = MTS_o_lo, file = "out_o/MTS_o_lo.csv")
# write.csv(x = MTS_o_hi, file = "out_o/MTS_o_hi.csv")
# write.csv(x = ddt.MTS_o, file = "out_o/ddt.MTS_o.csv")
# write.csv(x = ddt.MTS_o_lo, file = "out_o/ddt.MTS_o_lo.csv")
# write.csv(x = ddt.MTS_o_hi, file = "out_o/ddt.MTS_o_hi.csv")

#
###

###########################
## VISUALISE TIME SERIES ##
###########################

## Dark mode
par(bg="black", col.axis="white", col.lab="white", col.main="white", col="white")

## Plot MTS_o
num_series = length(unique(MTS[,1]))
num_variables = ncol(MTS)-2
par(mfrow=c(2,1), mar=c(3,3,1,1))
# par(mfrow=c(num_series,2), mar=c(3,3,1,1))
#
k = 2
for (k in 1:num_series)
{
  ## Time series index 
  idx = unique(MTS[,1])[k]
  
  ## Select time series
  s = (MTS[,1] == idx)
  TS = MTS[s,][,-1]
  TS_o = MTS_o[s,][,-1]
  TS_o_lo = MTS_o_lo[s,][,-1]
  TS_o_hi = MTS_o_hi[s,][,-1]
  ddt.TS_o = ddt.MTS_o[s,][,-1]
  ddt.TS_o_lo = ddt.MTS_o_lo[s,][,-1]
  ddt.TS_o_hi = ddt.MTS_o_hi[s,][,-1]
  
  ## Plot states
  colvect = rainbow(ncol(TS))
  for (i in 2:ncol(TS))
  {
    if (i == 2) plot(TS[,1], TS[,i], cex=0, ylim=c(-1,1)*3)
    lines(c(min(TS[,1]), max(TS[,1])), c(0,0), lty=2)
    x = TS_o[,1]
    y = TS_o[,i]
    y_lo = TS_o_lo[,i]
    y_hi = TS_o_hi[,i]
    polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[i])
    points(TS[,1], TS[,i], col=colvect[i], pch=16)
  }
  #
  ## Plot dynamics
  colvect = rainbow(ncol(TS))
  for (i in 2:ncol(TS))
  {
    if (i == 2) plot(TS[,1], TS[,i], cex=0, ylim=c(-1,1)*3)
    lines(c(min(TS[,1]), max(TS[,1])), c(0,0), lty=2)
    x = ddt.TS_o[,1]
    y = ddt.TS_o[,i]
    y_lo = ddt.TS_o_lo[,i]
    y_hi = ddt.TS_o_hi[,i]
    polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[i])
  }
}
#
par(mfrow=c(1,1))

#
###