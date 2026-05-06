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
source("f_bngm.r")
source("f_model_p.r")
source("f_utils.r")

#
###

#############
## MODEL 1 ##
#############

## Parameters of process model
N = 3 # 3 explanatory variables
W_p = 20
N_p = 2 * W_p * (2+N)
sd1_p = .1
sd2_p = c(1, 1)*1

## Indices
# s_v = 1
# s_r = 2:(N_p+1)
# s_p = (N_p+2):(2*N_p+1)

## functions
r = function(X, Omega) f_p.eval(X, Omega)
ddt.n = function(X, Omega) f_p.eval(X, Omega)
ddOmega.ddt.n = function(X, Omega) ddOmega.f_p.eval(X, Omega)
# ddz.r = function(X,Omega) t(ddx.f_p.eval(X, Omega))[,2]
# ddt.z = function(X,Omega) Omega[s_v]^2 * ddz.r(X, Omega[s_r]) + f_p.eval(X, Omega[s_p])

## Model
dLogPost = function(x, X__, Y__)
{
  logLik_ = -1/sd1_p^2 * sum((Y__[,1] - ddt.n(X__, x))^2)
  logPrior_ = -1/sd2_p[1]^2 * sum(x[s_r]^2)
  return(logLik_ + logPrior_)
}
ddOmega.dLogPost = function(x, X__, Y__)
{
  ddOmega.logLik_ = -1/sd1_p^2 * sum((Y__[,1] - ddt.n(X__, x))^2)
}

## Initiate parameters
Omega_0 = rnorm(N_p*2 + 1, 0, 0.001)

#
###