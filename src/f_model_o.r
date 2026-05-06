######
##  ##
######

## Goal:

## Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

#################################
## FUNCTIONS OBSERVATION MODEL ##
#################################

## f_o ##
## goal: compute predicted values of response variable at time step t
# t     - float  - time step 
# Omega - vector - parameters 
f_o = function(t,Omega)
{    
    Omega = matrix(Omega,ncol=3)
    return(t(Omega[,1])%*%sin(pi*(t*Omega[,2] + Omega[,3])))
}

## ddt.f_o ##
## goal: compute time derivative of the predicted response t time step t
# t     - float  - time step 
# Omega - vector - parameters 
ddt.f_o = function(t,Omega)
{    
    Omega = matrix(Omega,ncol=3)
    return(t(Omega[,1])%*%(pi*Omega[,2]*cos(pi*(t*Omega[,2] + Omega[,3]))))
}

## ddOmega.f_o ##
## goal: compute derivative of the predicted response wtr to each network parameter
# t     - float  - time step
# Omega - vector - parameters 
ddOmega.f_o = function(t,Omega)
{    
    Omega = matrix(Omega,ncol=3)
    dfdOmega_1 = sin(pi*(t*Omega[,2] + Omega[,3]))
    dfdOmega_2 = Omega[,1]*pi*t*cos(pi*(t*Omega[,2] + Omega[,3]))
    dfdOmega_3 = Omega[,1]*pi*1*cos(pi*(t*Omega[,2] + Omega[,3]))
    return(c(dfdOmega_1,dfdOmega_2,dfdOmega_3))
}

## *.eval ##
## goal: compute functions across multiple time steps
# t     - vector - time steps in arbitrary units
# Omega - vector - parameters 
f_o.eval = function(t,Omega) apply(t(t),2,function(x) f_o(x,Omega))
ddt.f_o.eval = function(t,Omega) apply(t(t),2,function(x) ddt.f_o(x,Omega))
ddOmega.f_o.eval = function(t,Omega) apply(t(t),2,function(x) ddOmega.f_o(x,Omega))

#
###

####################
## MAIN FUNCTIONS ##
####################

fit_model_o = function(TS, K_o=10, W_o=10, sd1_o=1.0, sd2_o=1.0, alpha_i=1)
{
  ## Variables
  N = ncol(TS) - 1
  n = nrow(TS)
  N_o = W_o*3
  
  ## Explanatory and response variable
  t = TS[,1]
  Y = TS[,-1]
  nt = seq(min(t), max(t), (t[2]-t[1])/alpha_i)
  
  ## Standardise time steps
  t_ = 1:length(t)
  nt_ = seq(min(t_), max(t_), (t_[2]-t_[1])/alpha_i)
  
  ## Standardise data
  Y_ = Y
  
  ## Interpolate time series
  Yhat_o = ddt.Yhat_o = NULL
  Yhat_o_lo = ddt.Yhat_o_lo = NULL
  Yhat_o_hi = ddt.Yhat_o_hi = NULL
  i = 1
  for(i in 1:N)
  {
    ## Anchor ensembling
    chain = NULL
    for(k in 1:K_o)
    {
      ## Fit
      Omega_0 = rnorm(N_o, 0, 0.001)
      # Omega_f = argmax.logPost(t_, Y_[,i], f_o.eval, ddOmega.f_o.eval, Omega_0, sd1_o, sd2_o)
      Omega_f = argmax.logPostLaplace(t_, Y_[,i], f_o.eval, ddOmega.f_o.eval, Omega_0, sd1_o, sd2_o)
      # Omega_f = argmax.logMarPost(t_, Y_[,i], f_o.eval, ddOmega.f_o.eval, Omega_0, 1.0)
      chain = rbind(chain, Omega_f)
    }
    
    ## Expectation
    Yhat_o_ = t(apply(chain, 1, FUN = function(x) f_o.eval(nt_, x)))
    Yhat_o_mean_ = apply(Yhat_o_, 2, mean)
    Yhat_o_lo_ = apply(Yhat_o_, 2, quantile, probs=0.1)
    Yhat_o_hi_ = apply(Yhat_o_, 2, quantile, probs=0.9)
    ddt.Yhat_o_ = t(apply(chain, 1, FUN = function(x) ddt.f_o.eval(nt_, x)))
    ddt.Yhat_o_mean_ = apply(ddt.Yhat_o_, 2, mean)
    ddt.Yhat_o_lo_ = apply(ddt.Yhat_o_, 2, quantile, probs=0.1)
    ddt.Yhat_o_hi_ = apply(ddt.Yhat_o_, 2, quantile, probs=0.9)
    
    ## Collect
    Yhat_o = cbind(Yhat_o, Yhat_o_mean_)
    Yhat_o_lo = cbind(Yhat_o_lo, Yhat_o_lo_)
    Yhat_o_hi = cbind(Yhat_o_hi, Yhat_o_hi_)
    ddt.Yhat_o = cbind(ddt.Yhat_o, ddt.Yhat_o_mean_)
    ddt.Yhat_o_lo = cbind(ddt.Yhat_o_lo, ddt.Yhat_o_lo_)
    ddt.Yhat_o_hi = cbind(ddt.Yhat_o_hi, ddt.Yhat_o_hi_)
  }
  
  ## Format
  Yhat_o = cbind(nt, Yhat_o)
  Yhat_o_lo = cbind(nt, Yhat_o_lo)
  Yhat_o_hi = cbind(nt, Yhat_o_hi)
  ddt.Yhat_o = cbind(nt, ddt.Yhat_o)
  ddt.Yhat_o_lo = cbind(nt, ddt.Yhat_o_lo)
  ddt.Yhat_o_hi = cbind(nt, ddt.Yhat_o_hi)
  
  ##
  colnames(Yhat_o) = colnames(TS)
  colnames(Yhat_o_lo) = colnames(TS)
  colnames(Yhat_o_hi) = colnames(TS)
  colnames(ddt.Yhat_o) = colnames(TS)
  colnames(ddt.Yhat_o_lo) = colnames(TS)
  colnames(ddt.Yhat_o_hi) = colnames(TS)
  
  ## Save
  return(list("Yhat_o"=Yhat_o, "Yhat_o_lo"=Yhat_o_lo, "Yhat_o_hi"=Yhat_o_hi, 
              "ddt.Yhat_o"=ddt.Yhat_o, "ddt.Yhat_o_lo"=ddt.Yhat_o_lo, "ddt.Yhat_o_hi"=ddt.Yhat_o_hi))

}

#
###