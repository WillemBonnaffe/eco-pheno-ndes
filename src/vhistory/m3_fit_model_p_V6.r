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

###############
## FUNCTIONS ##
###############

## Functions
ECI = function(X__, chain, func, lb=0.1, rb=0.9)
{
  Yhat_p__ = t(apply(chain, 1, FUN = function(x) func(X__, x)))
  Yhat_p_mean__ = apply(Yhat_p__, 2, median)
  Yhat_p_lo__ = apply(Yhat_p__, 2, quantile, probs=lb)
  Yhat_p_hi__ = apply(Yhat_p__, 2, quantile, probs=rb)
  return(list("mean" = Yhat_p_mean__, "lo" = Yhat_p_lo__, "hi" = Yhat_p_hi__))
}

#
###

##################
## FORMAT MTS_o ##
##################

##
MTS_o = data.frame(MTS_o)
ddt.MTS_o = data.frame(ddt.MTS_o)
for(i in 3:5) MTS_o[,i] = as.numeric(MTS_o[,i])
for(i in 3:5) ddt.MTS_o[,i] = as.numeric(ddt.MTS_o[,i])
head(MTS_o)
head(ddt.MTS_o)

#
###

#######################
## FIT PROCESS MODEL ##
#######################

## Train parameters
train_split = .75
train_lb = 0.25
train_rb = 0.75
t_lb = round(length(s_)*train_lb)
t_rb = round(length(s_)*train_rb)

## Split train and test
s_l = NULL
# selected_TS = c("D","J","K")
# selected_TS = c("A","D","F","G","J","K")
# selected_TS = c("B","E","H")
selected_TS = c("B","C","E","H","I","L")
# selected_TS = c("A","B","C","D","E","F","G","H","I","J","K","L")
s_bc = NULL
s_fc = NULL
for (selected_TS_ in selected_TS)
{
  ## Subset time series
  s_ = which((MTS_o[,1] == selected_TS_))
  
  ## Training set
  # s_l_ = s_[1:round(length(s_)*train_split)]
  s_l_ = s_[round(length(s_)*train_lb):round(length(s_)*train_rb)]
  
  ## Backcast and forecast set
  s_bc_ = s_[1:round(length(s_)*train_lb)]
  s_fc_ = s_[round(length(s_)*train_rb):length(s_)]
  
  ## Collect
  s_l = c(s_l, s_l_)
  s_bc = c(s_bc, s_bc_)
  s_fc = c(s_fc, s_fc_)
}
# s_l = which((MTS_o[,1] == "A") | MTS_o[,1] == "D")
# s_l = 1:round(length(s_l)*train_split)

## Variables
s = -c(1,2)
X = MTS_o[,s]
Y = ddt.MTS_o[,s]

## Standardise predictive variables
X_ = X
# mean_x = apply(X_[s_l,],2,mean)
# sd_x = apply(X_[s_l,],2,sd)
# X_ = t((t(X_)-mean_x)/sd_x)

## Standardise response variable
Y_ = Y
# mean_y = apply(Y_[s_l,],2,mean)
# sd_y = apply(Y_[s_l,],2,sd)
# Y_ = t((t(Y_))/sd_y) # not standardising wrt mean as 0 is informative

#
###

#############
## MODEL 1 ##
#############

## Parameters of process model
N = 3 # 3 explanatory variables
W_p = 10
N_p = 2 * W_p * (2+N)
sd1_p = .1
sd2_p = .7 # N_DJK = .25; Z_DJK = .45; N_BEH = 1.0; Z_BEH = 1.25
rv_idx = 1

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

#################
## PILOT CHAIN ##
#################

## Parameters
nIt = 3
sd2_p_vect = seq(0.1,1.0,0.1)

## Run chain
chain = NULL
chain_bc = NULL
chain_fc = NULL
for (sd2_p_ in sd2_p_vect)
{
  ## Iterator
  print(sd2_p_)
  
  for (k in 1:nIt)
  {
    ## Iterator
    print(k)
    
    ## Train
    Omega_0 = rnorm(N_p, 0, 0.001)
    chain_ = argmax.logPost(X=X_[s_l,],
                            Y=Y_[s_l,rv_idx],
                            f=ddt.n,
                            df=ddOmega.ddt.n,
                            Omega=Omega_0,
                            sd_1=sd1_p,
                            sd_2=sd2_p_)
    
    ## Test BC
    chain_bc_ = logLik(X=X_[s_bc,],
                       Y=Y_[s_bc,rv_idx],
                       f=ddt.n,
                       Omega=chain_,
                       sd_1=sd1_p
                       )

    ## Test FC
    chain_fc_ = logLik(X=X_[s_fc,],
                        Y=Y_[s_fc,rv_idx],
                        f=ddt.n,
                        Omega=chain_,
                        sd_1=sd1_p
                        )
    
    ## Collect
    chain = rbind(chain, chain_)
    chain_bc = rbind(chain_bc, chain_bc_)
    chain_fc = rbind(chain_fc, chain_fc_)
  }
}

## Format
chain_bc = matrix(chain_bc, ncol=nIt, byrow=T)
chain_fc = matrix(chain_fc, ncol=nIt, byrow=T)

#
###

###########################
## VISUALISE PILOT CHAIN ##
###########################

## Compute mean and quantiles
chain_bc_mean = apply(chain_bc, 1, mean)
chain_bc_q25 =  apply(chain_bc, 1, quantile, probs=0.25)
chain_bc_q75 =  apply(chain_bc, 1, quantile, probs=0.75)
#
chain_fc_mean = apply(chain_fc, 1, mean)
chain_fc_q25 =  apply(chain_fc, 1, quantile, probs=0.25)
chain_fc_q75 =  apply(chain_fc, 1, quantile, probs=0.75)

## Graphical parameters
par(mfrow=c(2,1))

## Backcast
x = sd2_p_vect
y = chain_bc_mean
y_lo = chain_bc_q25
y_hi = chain_bc_q75
plot(sd2_p_vect, chain_bc_mean, type="l", xlab="Regularisation Level (sd2_p)", ylab="Backcast Likelihood", ylim=c(min(y_lo), max(y_hi)))
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col="red")

## Forecast
x = sd2_p_vect
y = chain_fc_mean
y_lo = chain_fc_q25
y_hi = chain_fc_q75
plot(sd2_p_vect, chain_fc_mean, type="l", xlab="Regularisation Level (sd2_p)", ylab="Forecast Likelihood", ylim=c(min(y_lo), max(y_hi)))
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col="red")

## End
par(mfrow=c(1,1))

## Compromise
chain_bc_mean_std = (chain_bc_mean - min(chain_bc_mean))/(max(chain_bc_mean) - min(chain_bc_mean))
chain_fc_mean_std = (chain_fc_mean - min(chain_fc_mean))/(max(chain_fc_mean) - min(chain_fc_mean))
chain_bcfc_mean_std = (chain_fc_mean_std + chain_bc_mean_std)/2
x = sd2_p_vect
y = chain_bcfc_mean_std
plot(x, y, type="l", xlab="Regularisation Level (sd2_p)", ylab="Likelihood", ylim=c(0,1))
y = chain_bc_mean_std
lines(x, y, col="red", lty=2)
y = chain_fc_mean_std
lines(x, y, col="green", lty=2)
legend("bottom", legend = c("Backcast", "Forecast"), col=c("red", "green"), lty=2, horiz=T, bty="n")

#
###

###############
## RUN CHAIN ##
###############

## Paramaters
nIt = 30

## Run chain
chain = NULL
for (k in 1:nIt)
{
  print(k)
  Omega_0 = rnorm(N_p, 0, 0.001)
  chain_ = argmax.logPost(X=X_[s_l,],
                          Y=Y_[s_l,rv_idx],
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

for (selected_TS_ in selected_TS)
{
  ################################
  ## COMPUTE EXPECTATION AND CI ##
  ################################
  
  ## Get time series set
  s_t = which((MTS_o[,1] == selected_TS_))
  X__ = X_[s_t, ]
  Y__ = Y_[s_t, ]
  
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
  plot(x, rep(0,length(x)), ylim=c(-1,1)*4, lty=2, type="l", xlab="Time (by 6 Months)", ylab="Dynamics", main=paste("Time Series ", selected_TS_, sep=""))
  #
  y = ddt.N_p$mean
  y_lo = ddt.N_p$lo
  y_hi = ddt.N_p$hi
  polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
  lines(x, y, col=colvect[rv_idx])
  #
  ## Training data
  points(x, Y__[,rv_idx], col=colvect[rv_idx], pch=16)
  #
  ## Train bounds
  lines(c(t_lb, t_lb), c(-1,1)*4, lty=2)
  lines(c(t_rb, t_rb), c(-1,1)*4, lty=2)
  #
  ## Legend
  legend("bottom", legend = c("N", "Z", "E"), col=colvect, horiz=T, lty=1)
  
  ## Graphical parameters
  colvect = c("red", "green", "blue")
  
  ## Dynamics
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*4, lty=2, type="l", xlab="Time (by 6 Months)", ylab="Effects")
  #
  for (k in 1:N)
  {
    y = matrix(ddx.ddt.N_p$mean, ncol=3, byrow=T)[,k]
    y_lo = matrix(ddx.ddt.N_p$lo, ncol=3, byrow=T)[,k]
    y_hi = matrix(ddx.ddt.N_p$hi, ncol=3, byrow=T)[,k]
    polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[k])  
  }
  #
  ## Train bounds
  lines(c(t_lb, t_lb), c(-1,1)*4, lty=2)
  lines(c(t_rb, t_rb), c(-1,1)*4, lty=2)
  #
  ## Legend
  legend("bottom", legend = c("N", "Z", "E"), col=colvect, horiz=T, lty=1)
  
  ## End graph
  par(mfrow=c(1,1))
  
  #
  ###
}

################################
## COMPUTE EXPECTATION AND CI ##
################################

E = 0
for (E in c(-0.5, 0, 0.5))
{
  ## Graphical parameters
  par(bg="black", col.axis="white", col.lab="white", col.main="white", col="white")
  par(mfrow=c(3,1))
  colvect = c("red", "green", "blue")
  varnames = c("N", "Z", "E")
  
  for (s in 1:3)
  {
    ## Selected variable
    # s = 3
    
    ## Get x 
    x = seq(-3, 3, 0.1)# sort(rnorm(100,0,1))
    X__ = matrix(0,ncol=3,nrow=length(x))
    X__[,3] = E
    X__[,s] = x
    
    ## Expectation and credible interval
    ddt.N_p = ECI(X__, chain, ddt.n)
    ddt.N_p_interq = ECI(X__, chain, ddt.n, lb=0.25, rb=0.75)
    # ddt.Z_p = ECI(X__, chain, ddt.z)
    # ddt.Z_p_s = ECI(X__, chain, ddt.z_s)
    # ddt.Z_p_p = ECI(X__, chain, ddt.z_p)
    
    ## Dynamics
    x = X__[,s]
    plot(x, rep(0,length(x)), ylim=c(-1,1)*2, lty=2, type="l", xlab=varnames[s], ylab="Dynamics")
    #
    ## Credible interval
    y = ddt.N_p$mean
    y_lo = ddt.N_p$lo
    y_hi = ddt.N_p$hi
    polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[s])
    #
    ## Inter-quartile
    y = ddt.N_p_interq$mean
    y_lo = ddt.N_p_interq$lo
    y_hi = ddt.N_p_interq$hi
    polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[s])
    #
    ## Data
    for (selected_TS_ in selected_TS)
    {
      s_ts = which((MTS_o[,1] == selected_TS_))
      X__real = X_[s_ts, ]
      Y__real = Y_[s_ts, ]
      points(X__real[,s], Y__real[,rv_idx])  
    }
    #
    ## Legend
    legend("bottom", legend = varnames, col=colvect, horiz=T, lty=1, bty="n")
    legend("topright", legend = paste("E = ", E), bty="n")
  }
  #
  par(mfrow=c(1,1))
}

#
###