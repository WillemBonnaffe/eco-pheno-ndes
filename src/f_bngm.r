######
##  ##
######

## Goal:

## Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

##############
## INITIATE ##
##############

## Import modules
source("f_slp.r")

#
###

#####################################
## FUNCTIONS NORMAL BAYESIAN MODEL ##
#####################################

## goal: functions to define a simple Bayesian model with Gaussian error structure

## logLik ##
## goal: compute log likelihood of the process model 
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
# sd_1  - float  - standard deviation of likelihood
logLik = function(X,Y,f,Omega,sd_1)
{
  res = Y - f(X,Omega)
  logLik = - sum((res^2))/(sd_1^2)
  return(logLik)
}

## logPrior ##
## goal: compute log prior density of the process model 
# Omega - vector - parameters
# sd_2  - float  - standard deviation of prior
logPrior = function(Omega,sd_2)
{
  logPrior = - sum((Omega^2))/(sd_2^2)
  return(logPrior)
}

## logPost ##
## goal: log posterior distribution with normal error 
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# Omega - vector - parameters
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
logPost = function(X,Y,f,Omega,sd_1,sd_2)
{
    res = Y - f(X,Omega)
    logLik = - sum((res^2))/(sd_1^2)
    logPrior = - sum((Omega^2))/(sd_2^2)
    logPost = logLik + logPrior
    return(logPost)
}

## ddOmega.logPost ##
## goal: compute the derivate of the log posterior density wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
ddOmega.logPost = function(X,Y,f,df,Omega,sd_1,sd_2)
{
    res = Y - f(X,Omega)
    ddOmega.res = - df(X,Omega)
    ddOmega.logLik = - 2 * ddOmega.res%*%res/(sd_1^2)
    ddOmega.logPrior = - 2 * Omega/(sd_2^2)
    ddOmega.logPost = ddOmega.logLik + ddOmega.logPrior
    return(ddOmega.logPost)
}

## argmax.logPost ##
## goal: compute parameter vector that maximises log posterior density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - initial parameters 
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
argmax.logPost = function(X,Y,f,df,Omega,sd_1,sd_2)
{
    error_ = function(x) -logPost(X,Y,f,x,sd_1,sd_2)
    graderror_ = function(x) -ddOmega.logPost(X,Y,f,df,x,sd_1,sd_2)
    Omega = optim(par = Omega,
                  fn = error_,
                  gr = graderror_,
                  method = "BFGS"#,
                  # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
                  )$par
    return(Omega)
}

#
###

######################################
## FUNCTIONS LAPLACE BAYESIAN MODEL ##
######################################

## goal: functions to define a simple Bayesian model with Gaussian error structure

## logLikLaplace ##
## goal: compute log likelihood of the process model 
# X     - matrix - explanatory variables
# Y     - vector - response variable
# Omega - vector - parameters
# sd_1  - float  - standard deviation of likelihood
logLikLaplace = function(X,Y,f,Omega,sd_1)
{
  res = Y - f(X,Omega)
  logLik = - sum((res^2))/(sd_1^2)
  return(logLik)
}

## logPriorLaplace ##
## goal: compute log prior density of the process model 
# Omega - vector - parameters
# sd_2  - float  - standard deviation of prior
logPriorLaplace = function(Omega,sd_2)
{
  logPrior = - sum(abs(Omega))/(sd_2)
  return(logPrior)
}

## logPostLaplace ##
## goal: log posterior distribution with normal error 
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# Omega - vector - parameters
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
logPostLaplace = function(X,Y,f,Omega,sd_1,sd_2)
{
  res = Y - f(X,Omega)
  logLik = - sum((res^2))/(sd_1^2)
  logPrior = - sum(abs(Omega))/(sd_2)
  logPost = logLik + logPrior
  return(logPost)
}

## ddOmega.logPostLaplace ##
## goal: compute the derivate of the log posterior density wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
ddOmega.logPostLaplace = function(X,Y,f,df,Omega,sd_1,sd_2)
{
  res = Y - f(X,Omega)
  ddOmega.res = - df(X,Omega)
  ddOmega.logLik = - 2 * ddOmega.res%*%res/(sd_1^2)
  ddOmega.logPrior = - sign(Omega)/sd_2
  ddOmega.logPost = ddOmega.logLik + ddOmega.logPrior
  return(ddOmega.logPost)
}

## argmax.logPostLaplace ##
## goal: compute parameter vector that maximises log posterior density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - initial parameters 
# sd_1  - float  - standard deviation in likelihood 
# sd_2  - float  - standard deviation in priors 
argmax.logPostLaplace = function(X,Y,f,df,Omega,sd_1,sd_2)
{
  error_ = function(x) -logPostLaplace(X,Y,f,x,sd_1,sd_2)
  graderror_ = function(x) -ddOmega.logPostLaplace(X,Y,f,df,x,sd_1,sd_2)
  Omega = optim(par = Omega,
                fn = error_,
                gr = graderror_,
                method = "BFGS"#,
                # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
  )$par
  return(Omega)
}

#
###

#######################################
## FUNCTIONS MARGINAL BAYESIAN MODEL ##
#######################################

## logMarLik ##
## goal: compute the log marginal likelihood 
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
logMarLik = function(X,Y,f,Omega)
{
    res = Y - f(X,Omega)
    logMarLik = - 0.5 * length(Y) * log(0.5 * sum(res^2)   + 1)
    return(logMarLik)
}

## logMarPri ##
## goal: compute the log marginal prior density
# Omega - vector - parameters
logMarPri = function(Omega)
{
    logMarPri = - 0.5 * length(Omega) * log(0.5 * sum(Omega^2) + 1) 
    return(logMarPri)
}

## logMarPost ##
## goal: compute the log marginal posterior density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
# c     - scalar - parameter to control strength of regularisation
logMarPost = function(X,Y,f,Omega,c=1)
{
    res = Y - f(X,Omega)
    logMarLik = - 0.5 * length(Y) * log(0.5 * sum(res^2) + 1)
    logMarPri = - 0.5 * length(Omega) * log(0.5 * sum(Omega^2) + 1)
       logMarPos = logMarLik + c*logMarPri
    return(logMarPos)
}

## ddOmega.logMarPost ##
## goal: compute derivate of log marginal posterior density wtr to each parameter
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - parameters
# c     - scalar - parameter to control strength of regularisation
ddOmega.logMarPost = function(X,Y,f,df,Omega,c=1)
{
    res = Y - f(X,Omega)
    ddOmega.res = - df(X,Omega)
    ddOmega.logMarLik = - 0.5 * length(Y) * 1/(0.5 * sum(res^2) + 1) * 0.5 * ddOmega.res%*%res
    ddOmega.logMarPri = - 0.5 * length(Omega) * 1/(0.5 * sum(Omega^2) + 1) * Omega
    ddOmega.logMarPos = ddOmega.logMarLik + c*ddOmega.logMarPri
    return(ddOmega.logMarPos)
}

## argmax.logMarPost ##
## goal: compute parameter vector that maximises the log marginal density
# X     - matrix - explanatory variables
# Y     - vector - response variable 
# f     - func   - function to predict response
# df    - func   - function to compute gradient of predictive function 
# Omega - vector - initial parameters 
# c     - scalar - parameter to control strength of regularisation
argmax.logMarPost = function(X,Y,f,df,Omega,c=1)
{
    error_ = function(x) -logMarPost(X,Y,f,x,c)
    graderror_ = function(x) -ddOmega.logMarPost(X,Y,f,df,x,c)
    Omega = optim(par = Omega,
                  fn = error_,
                  gr = graderror_,
                  method = "BFGS"# ,
                  # control=list("trace"=1,"REPORT"=1,"maxit"=100) # for debugging
                  )$par
    return(Omega)
}

#
###