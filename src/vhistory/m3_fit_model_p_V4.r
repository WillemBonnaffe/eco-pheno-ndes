##############
## INITIATE ##
##############

## Local imports
source("f_betteRplots.r")
source("f_bngm.r")
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
sd1_p = .1
sd2_p = c(1, 1)*10
train_split = 1.0

## Split train and test
s_l = which((MTS_o[,1] == "3"))
# s_l = which((MTS_o[,1] == "1")|(MTS_o[,1] == "2")|(MTS_o[,1] == "3")|(MTS_o[,1] == "4"))
# s_l = s_l[1:round(length(s_l)*train_split)]
# s_l = 1:round(nrow(X_)*train_split)
s_t = -s_l

## Variables
s = -c(1,2)
X = MTS_o[,s]
Y = ddt.MTS_o[,s]

## Standardise predictive variables
X_ = X
mean_x = apply(X_[s_l,],2,mean)
sd_x = apply(X_[s_l,],2,sd)
X_ = t((t(X_)-mean_x)/sd_x)
# mean_x = apply(X_,2,mean)
# sd_x = apply(X_,2,sd)
# X_ = t((t(X_)-mean_x)/sd_x)

## Standardise response variable
Y_ = Y
mean_y = apply(Y_[s_l,],2,mean)
sd_y = apply(Y_[s_l,],2,sd)
Y_ = t((t(Y_))/sd_y)
# mean_y = apply(Y_,2,mean)
# sd_y = apply(Y_,2,sd)
# Y_ = t((t(Y_)-mean_y)/sd_y)

## Indices Model 1
s_v = 1
s_r = 2:(N_p+1)
s_p = (N_p+2):(2*N_p+1)

## functions
r = function(X,Omega) f_p.eval(X, Omega)
ddz.r = function(X,Omega) t(ddx.f_p.eval(X, Omega))[,2]
ddt.n = function(X,Omega) f_p.eval(X, Omega[s_r])
ddt.z = function(X,Omega) Omega[s_v]^2 * ddz.r(X, Omega[s_r]) + f_p.eval(X, Omega[s_p])
ddt.z_s = function(X,Omega) Omega[s_v]^2 * ddz.r(X, Omega[s_r])
ddt.z_p = function(X,Omega) f_p.eval(X, Omega[s_p])

## Model 3
target = function(x) 
{
  X__ = X_[s_l, ]
  Y__ = Y_[s_l, ]
  logLik_1_ = -sum((Y__[,1] - ddt.n(X__, x))^2) 
  logLik_2_ = -sum((Y__[,2] - ddt.z(X__, x))^2)
  logPrior_1_ = -1/sd2_p[1]^2 * sum(x[s_r]^2)
  logPrior_2_ = -1/sd2_p[2]^2 * sum(x[s_p]^2)
  return(logLik_1_ + logLik_2_)
}

# ## Model 2
# target = function(x) 
# {
#   X__ = X_[s_l, ]
#   Y__ = Y_[s_l, ]
#   logLik_1_ = logLik(X__, Y__[,1], ddt.n, x, sd1_p)
#   logLik_2_ = logLik(X__, Y__[,2], ddt.z, x, sd1_p)
#   logPrior_1_ = logPriorLaplace(x[s_r], sd2_p[1])*0
#   logPrior_2_ = logPriorLaplace(x[s_p], sd2_p[2])*0
#   return(logLik_1_ + logLik_2_ + logPrior_1_ + logPrior_2_)
# }

# ## Indices Model 2
# s_v = 1
# s_r = 2:(N_p+1)
# s_p = (N_p+2):(2*N_p+1)
# 
# ## functions
# r = function(X,Omega) f_p.eval(X, Omega)
# ddz.r = function(X,Omega) t(ddx.f_p.eval(X, Omega))[,2]
# ddt.n = function(X,Omega) f_p.eval(X, Omega[s_r])
# ddt.z = function(X,Omega) Omega[s_v]^2 * ddz.r(X, Omega[s_r]) + f_p.eval(X, Omega[s_p])
# ddt.z_s = function(X,Omega) Omega[s_v]^2 * ddz.r(X, Omega[s_r])
# ddt.z_p = function(X,Omega) f_p.eval(X, Omega[s_p])
# 
# ## Model 2
# target = function(x) 
# {
#   X__ = X_[s_l, ]
#   Y__ = Y_[s_l, ]
#   logLik_1_ = logLik(X__, Y__[,1], ddt.n, x, sd1_p)
#   logLik_2_ = logLik(X__, Y__[,2], ddt.z, x, sd1_p)
#   logPrior_1_ = logPriorLaplace(x[s_r], sd2_p[1])
#   logPrior_2_ = logPriorLaplace(x[s_p], sd2_p[2])
#   return(logLik_1_ + logLik_2_ + logPrior_1_ + logPrior_2_)
# }

## Run chain
nIt = 10000
Omega_0 = rnorm(N_p*2 + 1, 0, 0.001)
# Omega_0 = Omega_map
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

## Get time series set
X__ = X_[, ]
Y__ = Y_[, ]

## Contribution of selection and plasticity
ddt.N_p = ECI(X__, chain, ddt.n)
ddt.Z_p = ECI(X__, chain, ddt.z)
ddt.Z_p_s = ECI(X__, chain, ddt.z_s)
ddt.Z_p_p = ECI(X__, chain, ddt.z_p)

#
###

################################
## VISUALISE TIME SERIES - V2 ##
################################

## Graphical parameters
par(bg="black", col.axis="white", col.lab="white", col.main="white", col="white")
colvect = c("red", "blue")
par(mfrow=c(2,1))

## Select time series
s = MTS_o[,1] == 4

## N and Z dynamics
x = MTS_o[s,2]
plot(x, rep(0,length(x)), ylim=c(-1,1)*3, lty=2, type="l")
#
y = ddt.N_p$mean[s]
y_lo = ddt.N_p$lo[s]
y_hi = ddt.N_p$hi[s]
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col=colvect[1])
#
## Training data
points(x, Y__[s,1], col=colvect[1], pch=16)
#
y = ddt.Z_p$mean[s]
y_lo = ddt.Z_p$lo[s]
y_hi = ddt.Z_p$hi[s]
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col=colvect[2])
#
## Training data
points(x, Y__[s,2], col=colvect[2], pch=16)

## Z dynamics
x = MTS_o[s,2]
plot(x, rep(0,length(x)), ylim=c(-1,1)*3, lty=2, type="l")
colvect = c("magenta", "cyan")
#
y = ddt.Z_p_s$mean[s]
y_lo = ddt.Z_p_s$lo[s]
y_hi = ddt.Z_p_s$hi[s]
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col=colvect[1])
#
y = ddt.Z_p_p$mean[s]
y_lo = ddt.Z_p_p$lo[s]
y_hi = ddt.Z_p_p$hi[s]
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col=colvect[2])
#
## Training data
points(x, Y__[s,2])

## End graph
par(mfrow=c(1,1))

#
###