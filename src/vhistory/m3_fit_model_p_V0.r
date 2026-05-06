##############
## INITIATE ##
##############

source("f_betteRplots.r")
source("f_bngm.r")
source("f_model_o.r")
source("f_model_p.r")

#
###

##############
## INITIATE ##
##############

## Goal: load data, functions

## Load data
TS = read.table("data/TS_2.csv",sep=";",header=T)

## Extract time steps and columns of interest
selected_time_steps = 1:nrow(TS)
selected_columns = c("t","R","G","B")
TS = TS[selected_time_steps,]
TS = TS[,selected_columns]

## Shorten column names
column_names = c("times","A","Z","R")
colnames(TS) = column_names

## Transform variables
TS[,2] = log(TS[,2])
TS[,3] = log(TS[,3])
TS[,4] = log(TS[,4])

## Normalise time series
TS[,-1] = apply(TS[,-1], 2, function(x) (x-mean(x))/sd(x))

## Load results observation model
load("obj_results_o.RData")

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## Parameters of process model
K_p = 3
W_p = 10
N_p = 2 * W_p * (2+N)
sd1_p = 1.0
sd2_p = c(1, 1, 1)*0.25
train_split = 0.75

## Data specs
N = ncol(Yhat_o)

## Variables
X = Yhat_o
ddt.X = ddt.Yhat_o
Y = ddt.Yhat_o

## Standardise predictive variables
X_ = X
mean_x = apply(X_,2,mean)
sd_x = apply(X_,2,sd)
X_ = t((t(X_)-mean_x)/sd_x)

## Standardise response variable
Y_ = Y
mean_y = apply(Y_,2,mean)
sd_y = apply(Y_,2,sd)
Y_ = t((t(Y_)-mean_y)/sd_y)

## Split train and test
s_l = 1:round(nrow(X_)*train_split)
s_t = -s_l

## functions
Yhat = function(X,Omega) f_p.eval(X, Omega)
ddOmega.Yhat = function(X,Omega) ddOmega.f_p.eval(X, Omega)

## fit
Yhat_p = ddx.Yhat_p = Geber_p = NULL
Yhat_p_lo = ddx.Yhat_p_lo = Geber_p_lo = NULL
Yhat_p_hi = ddx.Yhat_p_hi = Geber_p_hi = NULL
for(i in 1:N)
{
    ## Anchor ensembling
    chain = NULL
    for(k in 1:K_p)
    {
      ## Fit
      Omega_0 = rnorm(N_p, 0, 0.001)
      Omega_f = argmax.logPost(X_[s_l,], Y_[s_l, i], Yhat, ddOmega.Yhat, Omega_0, sd1_p, sd2_p[i])
      chain = rbind(chain, Omega_f)
    }

    # ## Initiate chain
    # target = function(x) return(logPost(X_[s_l,], Y_[s_l, i], f_p.eval, x, sd1_p, sd2_p[i]))
    # nIt = 10000
    # Omega_0 = rnorm(N_p, 0, 0.001)
    # Omega_0 = argmax.logPost(X_[s_l,], Y_[s_l, i], Yhat, ddOmega.Yhat, Omega_0, sd1_p, sd2_p[i])
    # 
    # ## Run chain
    # chain = DEMCpp(list("dTarget"=target,
    #                     "Theta_0"=Omega_0,
    #                     "epsilon"=0.001,
    #                     "nIt"=nIt))$chainList
    # chain = chain[round(0.5*nIt):nIt, -1] # Burn
    # chain = chain[seq(1, nrow(chain), nrow(chain)/100), ] # Thin

    ## Expectation dynamics
    Yhat_p_ = t(apply(chain, 1, FUN = function(x) f_p.eval(X_, x)))
    Yhat_p_mean_ = apply(Yhat_p_, 2, mean)
    Yhat_p_lo_ = apply(Yhat_p_, 2, quantile, probs=0.1)
    Yhat_p_hi_ = apply(Yhat_p_, 2, quantile, probs=0.9)

    ## Expectation effects
    ddx.Yhat_p_ = t(apply(chain, 1, FUN = function(x) t(ddx.f_p.eval(X_, x))))
    ddx.Yhat_p_mean_ = t(apply(ddx.Yhat_p_, 2, mean))
    ddx.Yhat_p_lo_ = t(apply(ddx.Yhat_p_, 2, quantile, probs=0.1))
    ddx.Yhat_p_hi_ = t(apply(ddx.Yhat_p_, 2, quantile, probs=0.9))
    
    ## Expectation Geber
    Geber_p_ = t(apply(chain, 1, FUN = function(x) t(t(ddt.X) * ddx.f_p.eval(X_, x))))
    Geber_p_mean_ = t(apply(Geber_p_, 2, mean))
    Geber_p_lo_ = t(apply(Geber_p_, 2, quantile, probs=0.1))
    Geber_p_hi_ = t(apply(Geber_p_, 2, quantile, probs=0.9))
    
    ## Collect
    Yhat_p = cbind(Yhat_p, Yhat_p_mean_)
    Yhat_p_lo = cbind(Yhat_p_lo, Yhat_p_lo_)
    Yhat_p_hi = cbind(Yhat_p_hi, Yhat_p_hi_)
    #
    ddx.Yhat_p = cbind(ddx.Yhat_p, t(ddx.Yhat_p_mean_))
    ddx.Yhat_p_lo = cbind(ddx.Yhat_p_lo, t(ddx.Yhat_p_lo_))
    ddx.Yhat_p_hi = cbind(ddx.Yhat_p_hi, t(ddx.Yhat_p_hi_))
    #
    Geber_p = cbind(Geber_p, t(Geber_p_mean_))
    Geber_p_lo = cbind(Geber_p_lo, t(Geber_p_lo_))
    Geber_p_hi = cbind(Geber_p_hi, t(Geber_p_hi_))
}

## Save
save(list = c("X_",
              "Y_",
              "Yhat_p", 
              "Yhat_p_lo", 
              "Yhat_p_hi", 
              "ddx.Yhat_p", 
              "ddx.Yhat_p_lo", 
              "ddx.Yhat_p_hi",
              "Geber_p", 
              "Geber_p_lo", 
              "Geber_p_hi"), 
     file = "obj_results_p.RData")

#
###