#####
## ##
#####

## To do:
## - Think about latent variables
## - Think about re-introducing gradient-based fitting
## - Do artificial experiment
## - Re-balance data vs parameters in regularisation (N_p/(N_d + N_p))

##############
## INITIATE ##
##############

## Local imports
source("f_betteRplots.r")
source("f_bngm.r")
source("f_model_p.r")
source("f_utils.r")

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## Train parameters
train_split = 1.0

## Split train and test
s_l = which((MTS_o[,1] == "4"))
# s_l = 1:round(length(s_l)*train_split)
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

## Standardise response variable
Y_ = Y
mean_y = apply(Y_[s_l,],2,mean)
sd_y = apply(Y_[s_l,],2,sd)
Y_ = t((t(Y_))/sd_y) # not standardising wrt mean as 0 is informative

#
###

#############
## MODEL 1 ##
#############

## Parameters of process model
N = 3 # 3 explanatory variables
W_p = 10
N_p = 2 * W_p * (2+N)
sd1_p = 1
sd2_p = 1.0

## Indices
# s_v = 1
# s_r = 2:(N_p+1)
# s_p = (N_p+2):(2*N_p+1)

## functions
r = function(X, Omega) f_p.eval(X, Omega)
ddt.n = function(X, Omega) f_p.eval(X, Omega)
ddOmega.ddt.n = function(X, Omega) ddOmega.f_p.eval(X, Omega)
ddx.ddt.n = function(X, Omega) ddx.f_p.eval(X, Omega)
# ddz.r = function(X,Omega) t(ddx.f_p.eval(X, Omega))[,2]
# ddt.z = function(X,Omega) Omega[s_v]^2 * ddz.r(X, Omega[s_r]) + f_p.eval(X, Omega[s_p])
   
#
###

###############
## RUN CHAIN ##
###############

## Run chain
nIt = 3
# Omega_0 = Omega_map
chain = NULL
for (k in 1:nIt)
{
  print(k)
  ## Initiate parameters
  Omega_0 = rnorm(N_p, 0, 0.001)
  chain_ = argmax.logPost(X=X_[s_l,],
                          Y=Y_[s_l,1],
                          f=ddt.n,
                          df=ddOmega.ddt.n,
                          Omega=Omega_0,
                          sd_1=sd1_p,
                          sd_2=sd2_p)
  chain = rbind(chain, chain_)
}
Omega_map = chain[which.max(chain[,1]),]
# chain = chain[round(0.5*nIt):nIt, ] # Burn
# chain = chain[seq(1, nrow(chain), nrow(chain)/100), ] # Thin

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
X__ = X_[s_l, ]
Y__ = Y_[s_l, ]

## Contribution of selection and plasticity
ddt.N_p = ECI(X__, chain, ddt.n)
ddx.ddt.N_p = ECI(X__, chain, ddx.ddt.n)
# ddt.Z_p = ECI(X__, chain, ddt.z)
# ddt.Z_p_s = ECI(X__, chain, ddt.z_s)
# ddt.Z_p_p = ECI(X__, chain, ddt.z_p)

#
###

################################
## VISUALISE TIME SERIES - V2 ##
################################

## Graphical parameters
par(bg="black", col.axis="white", col.lab="white", col.main="white", col="white")
par(mfrow=c(2,1))

## Dynamics
x = 1:nrow(X__)
plot(x, rep(0,length(x)), ylim=c(-1,1)*3, lty=2, type="l")
#
y = ddt.N_p$mean
y_lo = ddt.N_p$lo
y_hi = ddt.N_p$hi
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col=colvect[1])
#
## Training data
points(x, Y__[,1], col=colvect[1], pch=16)

## Graphical parameters
colvect = c("red", "green", "blue")

## Dynamics
x = 1:nrow(X__)
plot(x, rep(0,length(x)), ylim=c(-1,1)*3, lty=2, type="l")
#
for (k in 1:N)
{
  y = matrix(ddx.ddt.N_p$mean, ncol=3, byrow=T)[,k]
  y_lo = matrix(ddx.ddt.N_p$lo, ncol=3, byrow=T)[,k]
  y_hi = matrix(ddx.ddt.N_p$hi, ncol=3, byrow=T)[,k]
  polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
  lines(x, y, col=colvect[k])  
}

## End graph
par(mfrow=c(1,1))

#
###