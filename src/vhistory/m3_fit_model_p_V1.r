##############
## INITIATE ##
##############

## Local imports
source("f_betteRplots.r")
source("f_bngm.r")
source("f_model_o.r")
source("f_model_p.r")
source("f_utils.r")

## Rcpp
library("Rcpp")
sourceCpp("f_demc.cpp")

#
###

##############
## INITIATE ##
##############

## Goal: load data, functions

## Time series

## Load data

## Load results observation model

## Remove time step column

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## Data specs
N = ncol(MTS_o) - 2

## Parameters of process model
K_p = 3
W_p = 10
N_p = 2 * W_p * (2+N)
sd1_p = .01
sd2_p = c(1, 1, 1)*1
train_split = 1.0

## Variables
s = -c(1,2)
X = MTS_o[,s]
Y = ddt.MTS_o[,s]

## Standardise predictive variables
X_ = X
# mean_x = apply(X_,2,mean)
# sd_x = apply(X_,2,sd)
# X_ = t((t(X_)-mean_x)/sd_x)

## Standardise response variable
Y_ = Y
# mean_y = apply(Y_,2,mean)
# sd_y = apply(Y_,2,sd)
# Y_ = t((t(Y_)-mean_y)/sd_y)

## Split train and test
s_l = which((MTS[,1] == "4"))
# s_l = which((MTS[,1] == "2")|(MTS[,1] == "4"))
# s_l = 1:round(nrow(X_)*train_split)
s_t = -s_l

## functions
r = function(X,Omega) f_p.eval(X, Omega)
ddz.r = function(X,Omega) t(ddx.f_p.eval(X, Omega))[,2]
ddt.n = function(X,Omega) f_p.eval(X, Omega)
ddt.z = function(X,Omega) Omega[1]^2 * ddz.r(X, Omega[(1+1):(N_p+1)]) + f_p.eval(X, Omega[(N_p+2):(2*N_p+1)])

## Model
s_ = 2:(N_p+1)
target = function(x) 
{
  X__ = X_[s_l, ]
  Y__ = Y_[s_l, ]
  target_1 = logPostLaplace(X__, Y__[,1], ddt.n, x[s_], sd1_p, sd2_p[1])
  target_2 = logPostLaplace(X__, Y__[,2], ddt.z, x, sd1_p, sd2_p[1])
  return(target_1 + target_2)
}

# ## Model 1
# target = function(x) 
# {
#   X__ = X_[s_l, ]
#   Y__ = Y_[s_l, ]
#   target_1 = logPost(X__, Y__[,1], r, x, sd1_p, sd2_p[1])
#   # target_2 = 0
#   target_2 = logPost(X__, Y__[,2], ddz.r, x, sd1_p, sd2_p[1])
#   return(target_1 + target_2)
# }

## Run chain
nIt = 10000
# Omega_0 = rnorm(N_p*2 + 1, 0, 0.001)
Omega_0 = Omega_map
chain = DEMCpp(list("dTarget"=target,
                    "Theta_0"=Omega_0,
                    "epsilon"=0.001,
                    "nIt"=nIt))$chainList
Omega_map = chain[which.max(chain[,1]),-1]
chain = chain[round(0.5*nIt):nIt, -1] # Burn
chain = chain[seq(1, nrow(chain), nrow(chain)/100), ] # Thin
Omega_mean = apply(chain[,-1], 2, mean)

#
###

##############################
## EXPECTATION COMPUTATIONS ##
##############################

## Initiate
Yhat_p = NULL
Yhat_p_lo = NULL
Yhat_p_hi = NULL

## Get time series set
X__ = X_[, ]
Y__ = Y_[, ]

## 
Yhat_p_mean_ = NULL
Yhat_p_lo_ = NULL
Yhat_p_hi_ = NULL

## Expectation dynamics (growth rate)
Yhat_p__ = t(apply(chain, 1, FUN = function(x) ddt.n(X__, x[s_])))
Yhat_p_mean__ = apply(Yhat_p__, 2, mean)
Yhat_p_lo__ = apply(Yhat_p__, 2, quantile, probs=0.1)
Yhat_p_hi__ = apply(Yhat_p__, 2, quantile, probs=0.9)
Yhat_p_mean_ = cbind(Yhat_p_mean_, Yhat_p_mean__)
Yhat_p_lo_ = cbind(Yhat_p_lo_, Yhat_p_lo__)
Yhat_p_hi_ = cbind(Yhat_p_hi_, Yhat_p_hi__)

## Expectation dynamics (phenotype)
Yhat_p__ = t(apply(chain, 1, FUN = function(x) ddt.z(X__, x)))
Yhat_p_mean__ = apply(Yhat_p__, 2, mean)
Yhat_p_lo__ = apply(Yhat_p__, 2, quantile, probs=0.1)
Yhat_p_hi__ = apply(Yhat_p__, 2, quantile, probs=0.9)
Yhat_p_mean_ = cbind(Yhat_p_mean_, Yhat_p_mean__)
Yhat_p_lo_ = cbind(Yhat_p_lo_, Yhat_p_lo__)
Yhat_p_hi_ = cbind(Yhat_p_hi_, Yhat_p_hi__)

## Expectation dynamics (environment)
Yhat_p_dummy__ = rep(NA, length(Yhat_p_mean__))
Yhat_p_mean_ = cbind(Yhat_p_mean_, Yhat_p_dummy__)
Yhat_p_lo_ = cbind(Yhat_p_lo_, Yhat_p_dummy__)
Yhat_p_hi_ = cbind(Yhat_p_hi_, Yhat_p_dummy__)
  
## Collect
Yhat_p = cbind(Yhat_p, Yhat_p_mean_)
Yhat_p_lo = cbind(Yhat_p_lo, Yhat_p_lo_)
Yhat_p_hi = cbind(Yhat_p_hi, Yhat_p_hi_)

#
###

#########################
##  FORMAT PREDICTIONS ##
#########################

## Add time series id and time step
ddt.MTS_p = cbind(MTS_o[,1:2], Yhat_p)
ddt.MTS_p_lo = cbind(MTS_o[,1:2], Yhat_p_lo)
ddt.MTS_p_hi = cbind(MTS_o[,1:2], Yhat_p_hi)

## Format
colnames(ddt.MTS_p) = colnames(MTS_o)

#
###

##################################
## DECOMPOSE PHENOTYPE DYNAMICS ##
##################################

## Functions
ECI = function(X__, chain, func)
{
  Yhat_p__ = t(apply(chain, 1, FUN = function(x) func(X__, x)))
  Yhat_p_mean__ = apply(Yhat_p__, 2, mean)
  Yhat_p_lo__ = apply(Yhat_p__, 2, quantile, probs=0.1)
  Yhat_p_hi__ = apply(Yhat_p__, 2, quantile, probs=0.9)
  return(list("mean" = Yhat_p_mean__, "lo" = Yhat_p_lo__, "hi" = Yhat_p_hi__))
}

##
ddt.z_s = function(X,Omega) Omega[1]^2 * ddz.r(X, Omega[(1+1):(N_p+1)])
ddt.z_p = function(X,Omega) f_p.eval(X, Omega[(N_p+2):(2*N_p+1)])

##
ddt.Z_p_s = ECI(X__, chain, ddt.z_s)
ddt.Z_p_p = ECI(X__, chain, ddt.z_p)

# ## Plot MTS_o vs MTS_p
# num_series = length(unique(MTS[,1]))
# num_variables = ncol(MTS)-2
# par(mfrow=c(num_series,4))
# #
# k = 2
# for (k in 1:num_series)
# {
#   ## Time series index 
#   idx = unique(MTS[,1])[k]
#   
#   ## Select time series
#   s = (MTS[,1] == idx)
#   TS = MTS[s,][,-1]
#   TS_o = MTS_o[s,][,-1]
#   TS_o_lo = MTS_o_lo[s,][,-1]
#   TS_o_hi = MTS_o_hi[s,][,-1]
#   ddt.TS_o = ddt.MTS_o[s,][,-1]
#   ddt.TS_o_lo = ddt.MTS_o_lo[s,][,-1]
#   ddt.TS_o_hi = ddt.MTS_o_hi[s,][,-1]
#   ddt.TS_p = ddt.MTS_p[s,][,-1]
#   ddt.TS_p_lo = ddt.MTS_p_lo[s,][,-1]
#   ddt.TS_p_hi = ddt.MTS_p_hi[s,][,-1]
#   
#   ## Plot selection vs plasticity
#   colvect = c("salmon", "pink") # rainbow(ncol(TS))
#   plot(TS[,1], TS[,i], cex=0, ylim=c(-1,1)*1)
#   lines(c(min(TS[,1]), max(TS[,1])), c(0,0), lty=2)
#   x = ddt.TS_p[,1] 
#   #
#   y = ddt.Z_p_s$mean[s]
#   y_lo = ddt.Z_p_s$lo[s]
#   y_hi = ddt.Z_p_s$hi[s]
#   polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
#   lines(x, y, col=colvect[1])
#   y = ddt.Z_p_p$mean[s]
#   y_lo = ddt.Z_p_p$lo[s]
#   y_hi = ddt.Z_p_p$hi[s]
#   polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
#   lines(x, y, col=colvect[2])
# }
# #
# par(mfrow=c(1,1))

#
###

###########################
## VISUALISE TIME SERIES ##
###########################

## Dark mode
par(bg="black", col.axis="white", col.lab="white", col.main="white", col="white")

## Plot MTS_o vs MTS_p
num_series = length(unique(MTS[,1]))
num_variables = ncol(MTS)-2
par(mfrow=c(num_series,4))
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
  ddt.TS_p = ddt.MTS_p[s,][,-1]
  ddt.TS_p_lo = ddt.MTS_p_lo[s,][,-1]
  ddt.TS_p_hi = ddt.MTS_p_hi[s,][,-1]
  
  ## Plot states
  colvect = c("","green", "purple", "red")# rainbow(ncol(TS))
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
  ## Plot observed dynamics
  colvect = c("","green", "purple", "red") # rainbow(ncol(TS))
  for (i in 2:ncol(TS))
  {
    if (i == 2) plot(TS[,1], TS[,i], cex=0, ylim=c(-1,1)*1)
    lines(c(min(TS[,1]), max(TS[,1])), c(0,0), lty=2)
    x = ddt.TS_o[,1]
    y = ddt.TS_o[,i]
    y_lo = ddt.TS_o_lo[,i]
    y_hi = ddt.TS_o_hi[,i]
    polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[i])
  }
  #
  ## Plot predicted dynamics
  colvect = c("","green", "purple", "red") # rainbow(ncol(TS))
  for (i in 2:ncol(TS))
  {
    if (i == 2) plot(TS[,1], TS[,i], cex=0, ylim=c(-1,1)*1)
    lines(c(min(TS[,1]), max(TS[,1])), c(0,0), lty=2)
    x = ddt.TS_p[,1]
    y = ddt.TS_p[,i]
    y_lo = ddt.TS_p_lo[,i]
    y_hi = ddt.TS_p_hi[,i]
    polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[i])
  }
  #
  ## Plot selection vs plasticity
  colvect = c("green", "purple") # rainbow(ncol(TS))
  plot(TS[,1], TS[,i], cex=0, ylim=c(-1,1)*1)
  lines(c(min(TS[,1]), max(TS[,1])), c(0,0), lty=2)
  x = ddt.TS_p[,1] 
  #
  y = ddt.Z_p_s$mean[s]
  y_lo = ddt.Z_p_s$lo[s]
  y_hi = ddt.Z_p_s$hi[s]
  polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
  lines(x, y, col=colvect[1])
  y = ddt.Z_p_p$mean[s]
  y_lo = ddt.Z_p_p$lo[s]
  y_hi = ddt.Z_p_p$hi[s]
  polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
  lines(x, y, col=colvect[2])
}
#
par(mfrow=c(1,1))

#
###
