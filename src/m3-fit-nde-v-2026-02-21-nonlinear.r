#####
## ##
#####

## Goal: Fit process model to time series.
## Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

## Note: Next try with simple time-differences.

## Update log:
## 2025-09-02 - Added collection of table of effects for ∆N(t) and ∆Z(t).
## 2026-02-21 - Added de-transformation of variables for visualisations.

pdf(file = "figures-2026-03-03-nonlinear.pdf", width = 12, height = 12)

##############
## INITIATE ##
##############

source("f_betteRplots.r")
source("f_bngm.r")
source("f_model_o.r")
source("f_model_p.r")
source("f_utils.r")

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

## General graphical parameters
par(bty='l')

#
###

###########################
## FIT OBSERVATION MODEL ##
###########################

k = 1

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
  MTS_o_ = MTS[MTS[,1] == idx,][,-1]
  # list_o = fit_model_o(TS, sd1_o=.1, sd2_o=.1, K_o=30)

  ## Compute difference on log scale
  ddt.MTS_o_ = MTS_o_
  ddt.MTS_o_[,c(2,3)] = log(ddt.MTS_o_[,c(2,3)])
  ddt.MTS_o_ = apply(ddt.MTS_o_, 2, diff)
  MTS_o_ = MTS_o_[1:nrow(MTS_o_)-1,] # Remove last time step as used to compute difference
  ddt.MTS_o_[,1] = MTS_o_[,1] # Set time to original value
  head(MTS_o_)
  head(ddt.MTS_o_)
    
  # ## Compute difference on natural scale
  # ddt.MTS_o_ = apply(MTS_o_, 2, diff)
  # MTS_o_ = MTS_o_[1:nrow(MTS_o_)-1,] # Remove last time step as used to compute difference
  # ddt.MTS_o_[,1] = MTS_o_[,1]
  # ddt.MTS_o_[,2] = as.numeric(ddt.MTS_o_[,2] / MTS_o_[,2])
  # ddt.MTS_o_[,3] = as.numeric(ddt.MTS_o_[,3] / MTS_o_[,3])
  # head(MTS_o_)
  # head(ddt.MTS_o_)
  
  ## Collect objects
  MTS_o = rbind(MTS_o, cbind(idx, MTS_o_))
  ddt.MTS_o = rbind(ddt.MTS_o, cbind(idx, ddt.MTS_o_))
  
}

## Format
ddt.MTS_o = data.frame(ddt.MTS_o)
for (i in 2:ncol(ddt.MTS_o)) ddt.MTS_o[,i] = as.numeric(ddt.MTS_o[,i])

## Compute mean
s = -c(1,2)
means = apply(MTS_o[,s], 2, mean, na.rm=T)
stds = apply(MTS_o[,s], 2, sd, na.rm=T)

## Standardise
s = -c(1,2)
MTS_o[,s] = t((t(MTS_o[,s])-means)/stds)

## Compute mean y
s = -c(1,2)
means_y = apply(ddt.MTS_o[,s], 2, mean, na.rm=T)
stds_y = apply(ddt.MTS_o[,s], 2, sd, na.rm=T)

## Standardise y
s = -c(1,2)
ddt.MTS_o[,s] = t((t(ddt.MTS_o[,s]*1))/stds_y)

# ## Save results
# system("mkdir out_o/")
# write.csv(x = MTS_o, file = "out_o/MTS_o.csv")
# write.csv(x = ddt.MTS_o, file = "out_o/ddt.MTS_o.csv")

#
###

###########################
## VISUALISE TIME SERIES ##
###########################

## Plot MTS_o
num_series = length(unique(MTS_o[,1]))
num_variables = ncol(MTS_o)-2
colvect = c('','blue', 'red', 'orange')
par(mfrow=c(6,4), mar=c(4.5,4.5,1,1), cex.lab=1.5)
#
k = 2
for (k in 1:num_series)
{
  ## Time series index 
  idx = unique(MTS_o[,1])[k]
  #
  ## Select time series
  s = (MTS_o[,1] == idx)
  TS_o = MTS_o[s,][,-1]
  ddt.TS_o = ddt.MTS_o[s,][,-1]
  #
  ## Plot states
  for (i in 2:3)
  {
    if (i == 2) plot(TS_o[,1], TS_o[,i], cex=0, ylim=c(-1,1)*3, xlab='Time (years)', ylab='Value (S.U.)', bty='l', xaxt='n', yaxt='n') 
    if (i == 2) add_axes_and_grid_at(TS_o[,1], TS_o[,i], at_x=1:6, at_y=-3:3, label=c(''))
    # lines(c(min(TS_o[,1]), max(TS_o[,1])), c(0,0), lty=2)
    x = TS_o[,1]
    y = TS_o[,i]
    lines(x, y, col=colvect[i])
    points(TS_o[,1], TS_o[,i], col=colvect[i], pch=16)
    legend('bottom', legend=c('N(t)', 'Z(t)'), col=colvect[2:3], horiz=T, bty='n', lty=1)
  }
  #
  ## Plot dynamics
  for (i in 2:3)
  {
    if (i == 2) plot(TS_o[,1], TS_o[,i], cex=0, ylim=c(-1,1)*3, xlab='Time (years)', ylab='Value (S.U.)', bty='l', xaxt='n', yaxt='n')
    if (i == 2) add_axes_and_grid_at(TS_o[,1], TS_o[,i], at_x=1:6, at_y=-3:3, label=c(''))
    lines(c(min(TS_o[,1]), max(TS_o[,1])), c(0,0), lty=2)
    x = as.numeric(ddt.TS_o[,1])
    y = as.numeric(ddt.TS_o[,i])
    lines(x, y, col=colvect[i])
    points(x, y, col=colvect[i], pch=16)
    polygon(x=c(x,rev(x)), y=c(rep(0,length(y)),rev(y)), col=adjustcolor(colvect[i], alpha=0.2), border=NA)
    # for (j in 1:length(x)) lines(c(x[j], x[j]), y=c(0, y[j]), lty=1, col=colvect[i]) # Add differences as bars
    legend('bottom', legend=c(expression(Delta * N(t)), expression(Delta * Z(t))), col=colvect[2:3], horiz=T, bty='n', lty=1)
  }
}
#
par(mfrow=c(1,1))

#
###

##################
## FORMAT MTS_o ##
##################

## format MTS_o
MTS_o = data.frame(MTS_o)
ddt.MTS_o = data.frame(ddt.MTS_o)
for(i in 3:ncol(MTS_o)) MTS_o[,i] = as.numeric(MTS_o[,i])
for(i in 3:ncol(MTS_o)) ddt.MTS_o[,i] = as.numeric(ddt.MTS_o[,i])
head(MTS_o)
head(ddt.MTS_o)

## Add fishing status
TS_fished = c("A","D","F","G","J","K") # Fished time series
idx_fished = multigrep(x = MTS_o$idx, TS_fished)
fishing_status = rep(0, nrow(MTS_o))
fishing_status[idx_fished] = 1
MTS_o$fishing_status = fishing_status
head(MTS_o)
tail(MTS_o)

## Add time of year
# MTS_o$census = rep(c(rep(c(0,1),6),0),12)
# head(MTS_o)
# tail(MTS_o)

#
###

##############################
## FORMAT DATA FOR TRAINING ##
##############################

## Train parameters
train_lb = 0.25
train_rb = 0.75
t_lb = round(nrow(MTS_o)*train_lb)
t_rb = round(nrow(MTS_o)*train_rb)

## Split train and test
# selected_TS = c("A","D","F","G","J","K") # Fished time series
# selected_TS = c("B","C","E","H","I","L") # Not fished time series
selected_TS = c("A","B","C","D","E","F","G","H","I","J","K","L") # All
s_l = NULL
s_bc = NULL
s_fc = NULL
for (selected_TS_ in selected_TS)
{
  ## Subset time series
  s_ = which((MTS_o[,1] == selected_TS_))
  
  ## Training set
  s_l_ = s_[round(length(s_)*train_lb):round(length(s_)*train_rb)]
  
  ## Backcast and forecast set
  s_bc_ = s_[1:round(length(s_)*train_lb)]
  s_fc_ = s_[round(length(s_)*train_rb):length(s_)]
  
  ## Collect
  s_l = c(s_l, s_l_)
  s_bc = c(s_bc, s_bc_)
  s_fc = c(s_fc, s_fc_)
}

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
N = 5 # number explanatory variables
W_p = 10
N_p = 2 * W_p * (2 + N)
sd1_p = .1
sd2_p = .1 # N_DJK = .25; Z_DJK = .45; N_BEH = 1.0; Z_BEH = 1.25
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
nIt = 10
sd2_p_vect = seq(0.01,0.5,0.02) # seq(0.1,1.0,0.1)

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
par(mfrow=c(3,1))

## Backcast
x = sd2_p_vect
y = chain_bc_mean
y_lo = chain_bc_q25
y_hi = chain_bc_q75
plot(sd2_p_vect, chain_bc_mean, type="l", xlab="Regularisation Level (sd2_p)", ylab="Backcast Likelihood", ylim=c(min(y_lo), max(y_hi)))
add_grid()
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col="red")

## Forecast
x = sd2_p_vect
y = chain_fc_mean
y_lo = chain_fc_q25
y_hi = chain_fc_q75
plot(sd2_p_vect, chain_fc_mean, type="l", xlab="Regularisation Level (sd2_p)", ylab="Forecast Likelihood", ylim=c(min(y_lo), max(y_hi)))
add_grid()
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col="red")

## Compromise
chain_bc_mean_std = (chain_bc_mean - min(chain_bc_mean))/(max(chain_bc_mean) - min(chain_bc_mean))
chain_fc_mean_std = (chain_fc_mean - min(chain_fc_mean))/(max(chain_fc_mean) - min(chain_fc_mean))
chain_bcfc_mean_std = (chain_fc_mean_std + chain_bc_mean_std)/2
x = sd2_p_vect
y = chain_bcfc_mean_std
plot(x, y, type="l", xlab="Regularisation Level (sd2_p)", ylab="Likelihood", ylim=c(0,1), bty='l')
add_grid()
lines(x, y)
y = chain_bc_mean_std
lines(x, y, col="red", lty=2)
y = chain_fc_mean_std
lines(x, y, col="green", lty=2)
legend("bottom", legend = c("Backcast", "Forecast"), col=c("red", "green"), lty=2, horiz=T, bty="n")

#
###

###########################
## SELECT REGULARISATION ##
###########################

sd2_p = .3

#
##

###############
## RUN CHAIN ##
###############

## Parameters
nIt = 100

## Run chain
chain = NULL
for (k in 1:nIt)
{
  print(k)
  Omega_0 = rnorm(N_p, 0, 0.001)
  chain_ = argmax.logPost(X=X_,# X=X_[s_l,],
                          Y=Y_[,rv_idx],# Y=Y_[s_l,rv_idx],
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

#######################
## VISUALISE EFFECTS ##
#######################

## Graphical parameters
par(mfrow=c(6,4), bty='l', mar=c(4.5,4.5,1,1), cex.lab=1.5)
colvect = c('blue', 'red', 'orange', 'cyan', 'magenta')

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
  
  ## Dynamics
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*4, lty=2, type="l", xlab="Time (years)", ylab="Dynamics")#, main=paste("Time Series ", selected_TS_, sep=""))
  add_grid()
  lines(x, rep(0,length(x)), lty=2)
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
  legend("bottom", legend = c(expression(Delta * N(t))), col=colvect[rv_idx], horiz=T, lty=1, bty='n')
  
  ## Effects
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*2, lty=2, type="l", xlab="Time (years)", ylab="Effects")
  add_grid()
  lines(x, rep(0,length(x)), lty=2)
  #
  for (k in 1:N)
  {
    y = matrix(ddx.ddt.N_p$mean, ncol=N, byrow=T)[,k]
    y_lo = matrix(ddx.ddt.N_p$lo, ncol=N, byrow=T)[,k]
    y_hi = matrix(ddx.ddt.N_p$hi, ncol=N, byrow=T)[,k]
    polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[k])
  }
  #
  ## Train bounds
  lines(c(t_lb, t_lb), c(-1,1)*2, lty=2)
  lines(c(t_rb, t_rb), c(-1,1)*2, lty=2)
  #
  ## Legend
  legend("bottom", legend = c("N:N", "Z:N", "T1:N", "T2:N", "F:N"), col=colvect, horiz=T, bty='n', pch=1)
  # legend("bottom", legend = c("∂(∆N)/∂N", "∂(∆N)/∂Z", "∂(∆N)/∂E", "∂(∆N)/∂F"), col=colvect, horiz=T, lty=1, bty='n')
  
  #
  ###
}

## End graph
par(mfrow=c(1,1))

#
###

################################
## COMPUTE SUMMARY OF EFFECTS ##
################################

summary_table_effects_on_delta_N = rbind(
  apply(matrix(ddx.ddt.N_p$mean, ncol=N, byrow=T), 2, mean),
  apply(matrix(ddx.ddt.N_p$lo, ncol=N, byrow=T), 2, mean),
  apply(matrix(ddx.ddt.N_p$hi, ncol=N, byrow=T), 2, mean)
)
summary_table_effects_on_delta_N = summary_table_effects_on_delta_N * stds_y[1] # Multiply by the scale of the response variable to de-standardise dynamics
summary_table_effects_on_delta_N = round(summary_table_effects_on_delta_N, 3)
summary_table_effects_on_delta_N = cbind(c("mid","low","high"),summary_table_effects_on_delta_N)
colnames(summary_table_effects_on_delta_N) = c("estimate","N:N", "Z:N", "T1:N", "T2:N", "F:N")
write.table(x=summary_table_effects_on_delta_N, file="tab-effects-population-nonlinear.csv", sep=",", row.names=F)

#
###

#############################
## VISUALISE CONTRIBUTIONS ##
#############################

## Graphical parameters
par(mfrow=c(6,4), bty='l', mar=c(4.5,4.5,1,1), cex.lab=1.5)
colvect = c('blue', 'red', 'orange', 'cyan', 'magenta')

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
  
  ## Contributions (Geber)
  dynamics = cbind(Y__, 0)
  jacobian = matrix(ddx.ddt.N_p$mean, ncol=N, byrow=T)
  geber = dynamics * jacobian
  
  #
  ###
  
  ################################
  ## VISUALISE TIME SERIES - V2 ##
  ################################
  
  ## Dynamics
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*4, lty=2, type="l", xlab="Time (years)", ylab="Dynamics")#, main=paste("Time Series ", selected_TS_, sep=""))
  add_grid()
  lines(x, rep(0,length(x)), lty=2)
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
  legend("bottom", legend = c(expression(Delta * N(t))), col=colvect[rv_idx], horiz=T, lty=1, bty='n')
  
  # ## Effects
  # x = 1:nrow(X__)
  # plot(x, rep(0,length(x)), ylim=c(-1,1)*2, lty=2, type="l", xlab="Time (years)", ylab="Effects")
  # add_grid()
  # lines(x, rep(0,length(x)), lty=2)
  # #
  # for (k in 1:N)
  # {
  #   y = matrix(ddx.ddt.N_p$mean, ncol=N, byrow=T)[,k]
  #   y_lo = matrix(ddx.ddt.N_p$lo, ncol=N, byrow=T)[,k]
  #   y_hi = matrix(ddx.ddt.N_p$hi, ncol=N, byrow=T)[,k]
  #   polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
  #   lines(x, y, col=colvect[k])
  # }
  # #
  # ## Train bounds
  # lines(c(t_lb, t_lb), c(-1,1)*2, lty=2)
  # lines(c(t_rb, t_rb), c(-1,1)*2, lty=2)
  # #
  # ## Legend
  # legend("bottom", legend = c("∂(∆N)/∂N", "∂(∆N)/∂Z", "∂(∆N)/∂E", "∂(∆N)/∂F"), col=colvect, horiz=T, lty=1, bty='n')
  
  ## Contributions
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*2, lty=2, type="l", xlab="Time (years)", ylab="Effects")
  add_grid()
  lines(x, rep(0,length(x)), lty=2)
  #
  for (k in 1:N)
  {
    y = geber[,k]
    lines(x, y, col=colvect[k])
    points(x, y, col=colvect[k], pch=16)
    # polygon(x=c(x,rev(x)), y=c(rep(0,length(y)),rev(y)), col=adjustcolor(colvect[k], alpha=0.2), border=NA)
    # for (j in 1:length(x)) lines(c(x[j], x[j]), y=c(0, y[j]), lty=1, col=colvect[i]) # Add differences as bars
  }
  #
  ## Train bounds
  lines(c(t_lb, t_lb), c(-1,1)*2, lty=2)
  lines(c(t_rb, t_rb), c(-1,1)*2, lty=2)
  #
  ## Legend
  legend("bottom", legend = c("N:N", "Z:N", "T1:N", "T2:N", "F:N"), col=colvect, horiz=T, bty='n', pch=1)
  # legend("bottom", legend = c("N:N", "Z:N", "E:N", "F:N"), col=colvect, horiz=T, lty=1, bty='n')
  
  #
  ###
}

## End graph
par(mfrow=c(1,1))

#
###

############################################
## COMPUTE EXPECTATION AND CI - VERSION 1 ##
############################################

## Graphical parameters
par(bg="white", col.axis="black", col.lab="black", col.main="black", col="black")
par(mfrow=c(N-1,N-1), mar=c(4.5,4.5,1,1), cex.lab=1.5)
colvect = c("lightblue", "salmon")  # Updated to distinguish F_ conditions
varnames = c("Log Density (S.U.)", "Log Body Mass (S.U.)", "Summer Temperature (S.U.)", "Winter Temperature (S.U.)")
pchvect = c(1,2)

## Plot
for (index_explanatory_variable in c(3, 4))
{
  for (E_ in c(-1.0, 1.5))
  {
    for (s in 1:(N-1))
    {
      ## Get x 
      x = seq(0, 1, 0.01) * (max(X[,s]) - min(X[,s])) + min(X[,s])
      plot(x, rep(0,length(x)), ylim=c(-1,1)*2.5, lty=2, type="l", xlab=varnames[s], ylab="P.c. Population Growth")
      add_grid()
      lines(x, rep(0,length(x)), lty=2)
      
      ## Data points
      for (selected_TS_ in selected_TS)
      {
        s_ts = which((MTS_o[,1] == selected_TS_))
        X__real = X_[s_ts, ]
        Y__real = Y_[s_ts, ]
        points(X__real[,s], Y__real[,rv_idx], pch=pchvect[X__real[,ncol(X_)]+1], col=c('black','red')[X__real[,ncol(X_)]+1])
      }
      
      ## Loop over Fishing conditions (F_ = 0 and F_ = 1)
      for (f_index in 1:2) {
        F_ = f_index - 1  # Convert index to 0 or 1
        
        ## Form X__
        X__ = matrix(0, ncol=N, nrow=length(x))
        X__[,index_explanatory_variable] = E_
        X__[,ncol(X_)] = F_
        X__[,s] = x
        
        ## Expectation and credible interval
        ddt.N_p = ECI(X__, chain, ddt.n)
        ddt.N_p_interq = ECI(X__, chain, ddt.n, lb=0.25, rb=0.75)
        
        ## Credible interval
        y = ddt.N_p$mean
        y_lo = ddt.N_p$lo
        y_hi = ddt.N_p$hi
        polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[f_index], 0.25))
        lines(x, y, col=colvect[f_index])
        
        ## Inter-quartile range
        y = ddt.N_p_interq$mean
        y_lo = ddt.N_p_interq$lo
        y_hi = ddt.N_p_interq$hi
        polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[f_index], 0.5))
        lines(x, y, col=colvect[f_index])
      }
      
      ## Legends
      legend("top", legend = paste("Temperature =", E_), bty="n")
      legend("bottom", legend = c("Not Harvested", "Harvested"), horiz=T, lty=1, col=colvect, bty="n")
    }
  } 
}
par(mfrow=c(1,1))

#
###

############################################
## COMPUTE EXPECTATION AND CI - VERSION 2 ##
############################################

## Graphical parameters
par(bg="white", col.axis="black", col.lab="black", col.main="black", col="black")
par(mfrow=c(2,2), mar=c(4.5,4.5,1,1), cex.lab=1.5)
layout(mat = rbind(c(1,2,5),c(3,4,6),7:9))
colvect = c("black", "red")  # Updated to distinguish F_ conditions
varnames = c("Density (fish/tank)", "Body Mass (mg)", "Summer Temperature (°C)", "Winter Temperature (°C)")
pchvect = c(1,2)
index_treatment = 5

## Plot
for (index_explanatory_variable in c(3))
{
  for (E_ in c(0.0))
  {
    for (s in 1:(N-1))
    {
      ## Get x 
      x = seq(0, 1, 0.01) * (max(X[,s]) - min(X[,s])) + min(X[,s])
      x_natural = x * stds[s] + means[s]  # De-standardise variables
      plot(x_natural, 
           rep(0,length(x)), ylim=c(-1,1)*2.5, lty=2, type="l", xlab=varnames[s], ylab="Population Growth")
      add_grid()
      lines(x_natural, 
            rep(0,length(x)), lty=2)
      
      ## Data points
      for (selected_TS_ in selected_TS)
      {
        s_ts = which((MTS_o[,1] == selected_TS_))
        X__real = X_[s_ts, ]
        Y__real = Y_[s_ts, ]
        points(X__real[,s] * stds[s] + means[s], # De-standardise variables
               Y__real[,rv_idx], pch=pchvect[X__real[,index_treatment]+1], col=c('black','red')[X__real[,index_treatment]+1])
      }
      
      ## Loop over Fishing conditions (F_ = 0 and F_ = 1)
      for (f_index in 1:2) {
        F_ = f_index - 1  # Convert index to 0 or 1
        
        ## Form X__
        X__ = matrix(0, ncol=N, nrow=length(x))
        X__[,index_explanatory_variable] = E_
        X__[,ncol(X_)] = F_
        X__[,s] = x
        
        ## Expectation and credible interval
        ddt.N_p = ECI(X__, chain, ddt.n)
        ddt.N_p_interq = ECI(X__, chain, ddt.n, lb=0.25, rb=0.75)
        
        ## Credible interval
        y = ddt.N_p$mean
        y_lo = ddt.N_p$lo
        y_hi = ddt.N_p$hi
        polygon(x=c(x_natural,rev(x_natural)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[f_index], 0.25))
        lines(x_natural, y, col=colvect[f_index])
        
        ## Inter-quartile range
        y = ddt.N_p_interq$mean
        y_lo = ddt.N_p_interq$lo
        y_hi = ddt.N_p_interq$hi
        polygon(x=c(x_natural,rev(x_natural)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[f_index], 0.5))
        lines(x_natural, y, col=colvect[f_index])
      }
      
      ## Legends
      # legend("top", legend = paste("Temperature =", E_), bty="n")
      legend("bottom", legend = c("Not Harvested", "Harvested"), horiz=T, lty=1, col=colvect, bty="n")
    }
  }
}
par(mfrow=c(1,1))

#
###

#######################
## COMPUTE R-SQUARED ##
#######################

## Compute residuals
predictions = ECI(X_, chain, ddt.n)
predictions = as.numeric(predictions$mean)
residuals = (Y_[,rv_idx] - predictions)

## Compute r2
r2 = 1 - sd(residuals)/sd(Y_[,rv_idx])
print(r2)
write.table(x=r2, file="tab-r2-population-nonlinear.csv", sep=",", row.names=F)

## QQ plot
sd_residuals = sd(residuals)
random_samples = rnorm(length(residuals), 0, sd_residuals)
plot(sort(random_samples), sort(residuals))
lines(c(-1,1)*3,c(-1,1)*3,lty=2)

## Visualise residuals
par(mfrow=c(1,1))
hist(residuals)
par(mfrow=c(1,1))

#
###

#######################
## VISUALISE HEATMAP ##
#######################

## Get mesh
res = 100
X_min = apply(X_, 2, min) * 1.25
X_max = apply(X_, 2, max) * 1.25
X_mid = apply(X_, 2, mean)
x1 = seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/res)
x2 = seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/res)
X1 = matrix(data=NA, nrow=res, ncol=res)
X2 = matrix(data=NA, nrow=res, ncol=res)
for(i in 1:res) {
  for (j in 1:res) {
    X1[i,j] = x1[i]
    X2[i,j] = x2[j]
  }
}

## Figure
par(mfrow=c(3,3))
#
predictions_N = list()
k = 1
for (X_5 in c(0,1)){
  for (X_4 in c(X_min[4], X_mid[4], X_max[4])){
    for (X_3 in c(X_min[3], X_mid[3], X_max[3])){
      
      ## Form mesh
      mesh = cbind(as.vector(X1), as.vector(X2), X_3, X_4, X_5)
      
      ## Predictions
      predictions = ECI(mesh, chain, ddt.n)
      predictions = predictions$mean
      predictions = matrix(predictions, nrow=res, ncol=res)
      predictions_N[[k]] = predictions
      
      ## Visualise mesh
      breaks = seq(-5,5,1.0)
      image(predictions, xaxt='n', yaxt='n', xlab='Density (fish/tank)', ylab='Body Mass  (mg)', col=hcl.colors(length(breaks)-1, "PiYG", rev=FALSE, alpha = 0.75), breaks=breaks)
      contour(predictions, add=T)
      contour(predictions, add=T, levels = c(0), col = 'red', lwd=2)
      x_at = seq(0,1,0.2)
      x1_labels = round(seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/100)[x_at * 100 + 1], 2) * stds[1] + means[1]
      x2_labels = round(seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/100)[x_at * 100 + 1], 2) * stds[2] + means[2]
      axis(side=1, at=x_at, labels=round(x1_labels))
      axis(side=2, at=x_at, labels=round(x2_labels))
      
      ## Add points
      x1_points_x = (X_[,1] - X_min[1]) / (X_max[1] - X_min[1]) 
      x2_points_x = (X_[,2] - X_min[2]) / (X_max[2] - X_min[2]) 
      points(x1_points_x, x2_points_x, pch=X_[,5]+1)
      
      ## Iterator
      k = k + 1
      
    }
  }
}
#
par(mfrow=c(1,1))

#
###

#############
## MODEL 2 ##
#############

## Parameters of process model
N = 5 # number explanatory variables
W_p = 10
N_p = 2 * W_p * (2 + N)
sd1_p = .1
sd2_p = .075 # N_DJK = .25; Z_DJK = .45; N_BEH = 1.0; Z_BEH = 1.25
rv_idx = 2

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
nIt = 10
sd2_p_vect = seq(0.01,0.5,0.02) # seq(0.1,1.0,0.1)

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
par(mfrow=c(3,1))

## Backcast
x = sd2_p_vect
y = chain_bc_mean
y_lo = chain_bc_q25
y_hi = chain_bc_q75
plot(sd2_p_vect, chain_bc_mean, type="l", xlab="Regularisation Level (sd2_p)", ylab="Backcast Likelihood", ylim=c(min(y_lo), max(y_hi)))
add_grid()
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col="red")

## Forecast
x = sd2_p_vect
y = chain_fc_mean
y_lo = chain_fc_q25
y_hi = chain_fc_q75
plot(sd2_p_vect, chain_fc_mean, type="l", xlab="Regularisation Level (sd2_p)", ylab="Forecast Likelihood", ylim=c(min(y_lo), max(y_hi)))
add_grid()
polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
lines(x, y, col="red")

## Compromise
chain_bc_mean_std = (chain_bc_mean - min(chain_bc_mean))/(max(chain_bc_mean) - min(chain_bc_mean))
chain_fc_mean_std = (chain_fc_mean - min(chain_fc_mean))/(max(chain_fc_mean) - min(chain_fc_mean))
chain_bcfc_mean_std = (chain_fc_mean_std + chain_bc_mean_std)/2
x = sd2_p_vect
y = chain_bcfc_mean_std
plot(x, y, type="l", xlab="Regularisation Level (sd2_p)", ylab="Likelihood", ylim=c(0,1), bty='l')
add_grid()
lines(x, y)
y = chain_bc_mean_std
lines(x, y, col="red", lty=2)
y = chain_fc_mean_std
lines(x, y, col="green", lty=2)
legend("bottom", legend = c("Backcast", "Forecast"), col=c("red", "green"), lty=2, horiz=T, bty="n")

#
###

###########################
## SELECT REGULARISATION ##
###########################

sd2_p = .4

#
##

###############
## RUN CHAIN ##
###############

## Paramaters
nIt = 100

## Run chain
chain = NULL
for (k in 1:nIt)
{
  print(k)
  Omega_0 = rnorm(N_p, 0, 0.001)
  chain_ = argmax.logPost(X=X_,# X=X_[s_l,],
                          Y=Y_[,rv_idx],# Y=Y_[s_l,rv_idx],
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

#######################
## VISUALISE EFFECTS ##
#######################

## Graphical parameters
par(mfrow=c(6,4), bty='l', mar=c(4.5,4.5,1,1), cex.lab=1.5)
colvect = c('blue', 'red', 'orange', 'cyan', 'magenta')

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
  
  ## Dynamics
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*4, lty=2, type="l", xlab="Time (years)", ylab="Dynamics")#, main=paste("Time Series ", selected_TS_, sep=""))
  add_grid()
  lines(x, rep(0,length(x)), lty=2)
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
  legend("bottom", legend = c(expression(Delta * Z(t))), col=colvect[rv_idx], horiz=T, lty=1, bty='n')
  
  ## Effects
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*2, lty=2, type="l", xlab="Time (years)", ylab="Effects")
  add_grid()
  lines(x, rep(0,length(x)), lty=2)
  #
  for (k in 1:N)
  {
    y = matrix(ddx.ddt.N_p$mean, ncol=N, byrow=T)[,k]
    y_lo = matrix(ddx.ddt.N_p$lo, ncol=N, byrow=T)[,k]
    y_hi = matrix(ddx.ddt.N_p$hi, ncol=N, byrow=T)[,k]
    polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[k])
  }
  #
  ## Train bounds
  lines(c(t_lb, t_lb), c(-1,1)*2, lty=2)
  lines(c(t_rb, t_rb), c(-1,1)*2, lty=2)
  #
  ## Legend
  legend("bottom", legend = c("N:Z", "Z:Z", "T1:Z", "T2:Z", "F:Z"), col=colvect, horiz=T, bty='n', pch=1)
  # legend("bottom", legend = c("∂(∆N)/∂N", "∂(∆N)/∂Z", "∂(∆N)/∂E1", "∂(∆N)/∂E2", "∂(∆N)/∂F"), col=colvect, horiz=T, lty=1, bty='n')
  
  #
  ###
}

## End graph
par(mfrow=c(1,1))

#
###

################################
## COMPUTE SUMMARY OF EFFECTS ##
################################

summary_table_effects_on_delta_Z = rbind(
  apply(matrix(ddx.ddt.N_p$mean, ncol=N, byrow=T), 2, mean),
  apply(matrix(ddx.ddt.N_p$lo, ncol=N, byrow=T), 2, mean),
  apply(matrix(ddx.ddt.N_p$hi, ncol=N, byrow=T), 2, mean)
)
summary_table_effects_on_delta_Z = summary_table_effects_on_delta_Z * stds_y[2] # Multiply by the scale of the response variable 
summary_table_effects_on_delta_Z = round(summary_table_effects_on_delta_Z, 3)
summary_table_effects_on_delta_Z = cbind(c("mid","low","high"),summary_table_effects_on_delta_Z)
colnames(summary_table_effects_on_delta_Z) = c("estimate","N:Z", "Z:Z", "T1:Z", "T2:Z", "F:Z")
write.table(x=summary_table_effects_on_delta_Z, file="tab-effects-phenotype-nonlinear.csv", sep=",", row.names=F)

#
###

#############################
## VISUALISE CONTRIBUTIONS ##
#############################

## Graphical parameters
par(mfrow=c(6,4), bty='l', mar=c(4.5,4.5,1,1), cex.lab=1.5)
colvect = c('blue', 'red', 'orange','cyan', 'magenta')

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
  
  ## Contributions (Geber)
  dynamics = cbind(Y__, 0)
  jacobian = matrix(ddx.ddt.N_p$mean, ncol=N, byrow=T)
  geber = dynamics * jacobian
  
  #
  ###
  
  ################################
  ## VISUALISE TIME SERIES - V2 ##
  ################################
  
  ## Dynamics
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*4, lty=2, type="l", xlab="Time (years)", ylab="Dynamics")#, main=paste("Time Series ", selected_TS_, sep=""))
  add_grid()
  lines(x, rep(0,length(x)), lty=2)
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
  legend("bottom", legend = c(expression(Delta * Z(t))), col=colvect[rv_idx], horiz=T, lty=1, bty='n')
  
  # ## Effects
  # x = 1:nrow(X__)
  # plot(x, rep(0,length(x)), ylim=c(-1,1)*2, lty=2, type="l", xlab="Time (years)", ylab="Effects")
  # add_grid()
  # lines(x, rep(0,length(x)), lty=2)
  # #
  # for (k in 1:N)
  # {
  #   y = matrix(ddx.ddt.N_p$mean, ncol=N, byrow=T)[,k]
  #   y_lo = matrix(ddx.ddt.N_p$lo, ncol=N, byrow=T)[,k]
  #   y_hi = matrix(ddx.ddt.N_p$hi, ncol=N, byrow=T)[,k]
  #   polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
  #   lines(x, y, col=colvect[k])
  # }
  # #
  # ## Train bounds
  # lines(c(t_lb, t_lb), c(-1,1)*2, lty=2)
  # lines(c(t_rb, t_rb), c(-1,1)*2, lty=2)
  # #
  # ## Legend
  # legend("bottom", legend = c("∂(∆N)/∂N", "∂(∆N)/∂Z", "∂(∆N)/∂E", "∂(∆N)/∂F"), col=colvect, horiz=T, lty=1, bty='n')
  
  ## Contributions
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*2, lty=2, type="l", xlab="Time (years)", ylab="Effects")
  add_grid()
  lines(x, rep(0,length(x)), lty=2)
  #
  for (k in 1:N)
  {
    y = geber[,k]
    lines(x, y, col=colvect[k])
    points(x, y, col=colvect[k], pch=16)
    # polygon(x=c(x,rev(x)), y=c(rep(0,length(y)),rev(y)), col=adjustcolor(colvect[k], alpha=0.2), border=NA)
    # for (j in 1:length(x)) lines(c(x[j], x[j]), y=c(0, y[j]), lty=1, col=colvect[i]) # Add differences as bars
  }
  #
  ## Train bounds
  lines(c(t_lb, t_lb), c(-1,1)*2, lty=2)
  lines(c(t_rb, t_rb), c(-1,1)*2, lty=2)
  #
  ## Legend
  legend("bottom", legend = c("N:Z", "Z:Z", "T1:Z", "T2:Z", "F:Z"), col=colvect, horiz=T, bty='n', pch=1)
  
  #
  ###
}

## End graph
par(mfrow=c(1,1))

#
###

################################
## COMPUTE EXPECTATION AND CI ##
################################

## Graphical parameters
par(bg="white", col.axis="black", col.lab="black", col.main="black", col="black")
par(mfrow=c(N-1,N-1), mar=c(4.5,4.5,1,1), cex.lab=1.5)
colvect = c("lightblue", "salmon")  # Updated to distinguish F_ conditions
varnames = c("Log Density (S.U.)", "Log Body Mass (S.U.)", "Summer Temperature (S.U.)", "Winter Temperature")
pchvect = c(1,2)
index_treatment = 5

## Plot
for (index_explanatory_variable in c(3,4))
{
  for (E_ in c(-1.0, 1.5))
  {
    for (s in 1:(N-1))
    {
      ## Get x 
      x = seq(0, 1, 0.01) * (max(X[,s]) - min(X[,s])) + min(X[,s])
      plot(x, rep(0,length(x)), ylim=c(-1,1)*2.5, lty=2, type="l", xlab=varnames[s], ylab="Phenotypic Change")
      add_grid()
      lines(x, rep(0,length(x)), lty=2)
      
      ## Data points
      for (selected_TS_ in selected_TS)
      {
        s_ts = which((MTS_o[,1] == selected_TS_))
        X__real = X_[s_ts, ]
        Y__real = Y_[s_ts, ]
        points(X__real[,s], Y__real[,rv_idx], pch=pchvect[X__real[,index_treatment]+1], col=c('black','red')[X__real[,index_treatment]+1])
      }
      
      ## Loop over Fishing conditions (F_ = 0 and F_ = 1)
      for (f_index in 1:2) {
        F_ = f_index - 1  # Convert index to 0 or 1
        
        ## Form X__
        X__ = matrix(0, ncol=N, nrow=length(x))
        X__[,index_explanatory_variable] = E_
        X__[,ncol(X_)] = F_
        X__[,s] = x
        
        ## Expectation and credible interval
        ddt.N_p = ECI(X__, chain, ddt.n)
        ddt.N_p_interq = ECI(X__, chain, ddt.n, lb=0.25, rb=0.75)
        
        ## Credible interval
        y = ddt.N_p$mean
        y_lo = ddt.N_p$lo
        y_hi = ddt.N_p$hi
        polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[f_index], 0.25))
        lines(x, y, col=colvect[f_index])
        
        ## Inter-quartile range
        y = ddt.N_p_interq$mean
        y_lo = ddt.N_p_interq$lo
        y_hi = ddt.N_p_interq$hi
        polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[f_index], 0.5))
        lines(x, y, col=colvect[f_index])
      }
      
      ## Legends
      legend("top", legend = paste("Temperature =", E_), bty="n")
      legend("bottom", legend = c("Not Harvested", "Harvested"), horiz=T, lty=1, col=colvect, bty="n")
    }
  }
}
par(mfrow=c(1,1))

#
###

############################################
## COMPUTE EXPECTATION AND CI - VERSION 2 ##
############################################

## Graphical parameters
par(bg="white", col.axis="black", col.lab="black", col.main="black", col="black")
par(mfrow=c(2,2), mar=c(4.5,4.5,1,1), cex.lab=1.5)
layout(mat = rbind(c(1,2,5),c(3,4,6),7:9))
colvect = c("black", "red")  # Updated to distinguish F_ conditions
varnames = c("Density (fish/tank)", "Body Mass (mg)", "Summer Temperature (°C)", "Winter Temperature (°C)")
pchvect = c(1,2)
index_treatment = 5

## Plot
for (index_explanatory_variable in c(3))
{
  for (E_ in c(0.0))
  {
    for (s in 1:(N-1))
    {
      ## Get x 
      x = seq(0, 1, 0.01) * (max(X[,s]) - min(X[,s])) + min(X[,s])
      x_natural = x * stds[s] + means[s]  # De-standardise variables
      plot(x_natural, 
           rep(0,length(x)), ylim=c(-1,1)*2.5, lty=2, type="l", xlab=varnames[s], ylab="Phenotypic Change")
      add_grid()
      lines(x_natural, 
            rep(0,length(x)), lty=2)
      
      ## Data points
      for (selected_TS_ in selected_TS)
      {
        s_ts = which((MTS_o[,1] == selected_TS_))
        X__real = X_[s_ts, ]
        Y__real = Y_[s_ts, ]
        points(X__real[,s] * stds[s] + means[s], # De-standardise variables
        Y__real[,rv_idx], pch=pchvect[X__real[,index_treatment]+1], col=c('black','red')[X__real[,index_treatment]+1])
      }
      
      ## Loop over Fishing conditions (F_ = 0 and F_ = 1)
      for (f_index in 1:2) {
        F_ = f_index - 1  # Convert index to 0 or 1
        
        ## Form X__
        X__ = matrix(0, ncol=N, nrow=length(x))
        X__[,index_explanatory_variable] = E_
        X__[,ncol(X_)] = F_
        X__[,s] = x
        
        ## Expectation and credible interval
        ddt.N_p = ECI(X__, chain, ddt.n)
        ddt.N_p_interq = ECI(X__, chain, ddt.n, lb=0.25, rb=0.75)
        
        ## Credible interval
        y = ddt.N_p$mean
        y_lo = ddt.N_p$lo
        y_hi = ddt.N_p$hi
        polygon(x=c(x_natural,rev(x_natural)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[f_index], 0.25))
        lines(x_natural, y, col=colvect[f_index])
        
        ## Inter-quartile range
        y = ddt.N_p_interq$mean
        y_lo = ddt.N_p_interq$lo
        y_hi = ddt.N_p_interq$hi
        polygon(x=c(x_natural,rev(x_natural)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[f_index], 0.5))
        lines(x_natural, y, col=colvect[f_index])
      }
      
      ## Legends
      # legend("top", legend = paste("Temperature =", E_), bty="n")
      legend("bottom", legend = c("Not Harvested", "Harvested"), horiz=T, lty=1, col=colvect, bty="n")
    }
  }
}
par(mfrow=c(1,1))

#
###

#######################
## COMPUTE R-SQUARED ##
#######################

## Compute residuals
predictions = ECI(X_, chain, ddt.n)
predictions = as.numeric(predictions$mean)
residuals = (Y_[,rv_idx] - predictions)

## Compute r2
r2 = 1 - sd(residuals)/sd(Y_[,rv_idx])
print(r2)
write.table(x=r2, file="tab-r2-phenotype-nonlinear.csv", sep=",", row.names=F)

## QQ plot
sd_residuals = sd(residuals)
random_samples = rnorm(length(residuals), 0, sd_residuals)
plot(sort(random_samples), sort(residuals))
lines(c(-1,1)*3,c(-1,1)*3,lty=2)

## Visualise residuals
par(mfrow=c(1,1))
hist(residuals)
par(mfrow=c(1,1))

#
###

#######################
## VISUALISE HEATMAP ##
#######################

## Get mesh
res = 100
X_min = apply(X_, 2, min) * 1.25
X_max = apply(X_, 2, max) * 1.25
X_mid = apply(X_, 2, mean)
x1 = seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/res)
x2 = seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/res)
X1 = matrix(data=NA, nrow=res, ncol=res)
X2 = matrix(data=NA, nrow=res, ncol=res)
for(i in 1:res) {
  for (j in 1:res) {
    X1[i,j] = x1[i]
    X2[i,j] = x2[j]
  }
}

## Figure
par(mfrow=c(3,3))
#
predictions_Z = list()
k = 1
for (X_5 in c(0,1)){
  for (X_4 in c(X_min[4], X_mid[4], X_max[4])){
    for (X_3 in c(X_min[3], X_mid[3], X_max[3])){
      
      ## Form mesh
      mesh = cbind(as.vector(X1), as.vector(X2), X_3, X_4, X_5)
      
      ## Predictions
      predictions = ECI(mesh, chain, ddt.n)
      predictions = predictions$mean
      predictions = matrix(predictions, nrow=res, ncol=res)
      predictions_Z[[k]] = predictions
      
      ## Visualise mesh
      breaks = seq(-5,5,1.0)
      image(predictions, xaxt='n', yaxt='n', xlab='Density (fish/tank)', ylab='Body Mass (mg)', col=hcl.colors(length(breaks)-1, "PiYG", rev=FALSE, alpha = 0.75), breaks=breaks)
      contour(predictions, add=T)
      contour(predictions, add=T, levels = c(0), col = 'red', lwd=2)
      x_at = seq(0,1,0.2)
      x1_labels = round(seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/100)[x_at * 100 + 1], 2) * stds[1] + means[1]
      x2_labels = round(seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/100)[x_at * 100 + 1], 2) * stds[2] + means[2]
      axis(side=1, at=x_at, labels=round(x1_labels))
      axis(side=2, at=x_at, labels=round(x2_labels))
      
      ## Add points
      x1_points_x = (X_[,1] - X_min[1]) / (X_max[1] - X_min[1]) 
      x2_points_x = (X_[,2] - X_min[2]) / (X_max[2] - X_min[2]) 
      points(x1_points_x, x2_points_x, pch=X_[,5]+1)
      
      ## Iterator
      k = k + 1
      
    }
  }
}
#
par(mfrow=c(1,1))

#
###

#####################
## LOCAL FUNCTIONS ##
#####################

add_streamlines = function()
{
  for (point in start_points) {
    
    # Initialize position
    x_pos = point[1]
    y_pos = point[2]
    cum_dist = 0  # Track cumulative distance traveled
    arrow_spacing = 0.15  # Adjust this to control arrow spacing
    
    for (step in 1:num_steps) {
      # Find the closest grid index
      i = round(x_pos * (res - 1)) + 1
      j = round(y_pos * (res - 1)) + 1
      
      # Ensure indices are within valid range
      if (i < 1 || i > res || j < 1 || j > res) break
      
      # Get vector field values
      dx = predictions_N[[k]][i, j] * step_size
      dy = predictions_Z[[k]][i, j] * step_size
      
      # Compute movement distance
      step_dist = sqrt(dx^2 + dy^2)
      cum_dist = cum_dist + step_dist  # Accumulate movement distance
      
      # Compute new position
      x_new = x_pos + dx
      y_new = y_pos + dy
      
      # Ensure new position stays within bounds
      if (x_new < 0 || x_new > 1 || y_new < 0 || y_new > 1) break
      
      # Only place an arrow if cumulative distance exceeds threshold
      if (cum_dist >= arrow_spacing) {
        arrows(x_pos, y_pos, x_new, y_new, length=0.05, col=adjustcolor("gray20", alpha=1.0), lwd=1.0)
        cum_dist = 0  # Reset cumulative distance after placing an arrow
      } else {
        segments(x_pos, y_pos, x_new, y_new, col=adjustcolor("gray20", alpha=1.0), lwd=1.0)
      }
      
      # Update position
      x_pos = x_new
      y_pos = y_new
    }
  }
}

#
###

#######################
## COMBINED ANALYSIS ##
#######################

## Graphical parameters
breaks = c(1,2,3,4,5)

index = NULL
i = 1; j = 1:3; k = 2
index = c(index, (i-1)*9 + (j-1)*3 + k)
i = 2; j = 1:3; k = 2
index = c(index, (i-1)*9 + (j-1)*3 + k)

## Figure
par(mfrow=c(3,3), bty='o')
#
for(k in 1:length(predictions_N)){
  
  ## Plot
  img = (predictions_N[[k]] > 0) + (predictions_Z[[k]] < 0) * 2
  image(img, col=hcl.colors(length(breaks)-1, "RdYlBu", rev=FALSE, alpha = 0.3), xaxt='n', yaxt='n', xlab='Density (fish/tank)', ylab='Body Mass (mg)')
  # image(img, col=hcl.colors(length(breaks)-1, "RdBu", rev=FALSE, alpha = 0.5), xaxt='n', yaxt='n', xlab='Log Density', ylab='Log Body Mass')
  # contour((predictions_N[[k]]) * (predictions_Z[[k]]), levels=0, add=T, col='red', lwd=4)  
  
  ## Axes
  x_at = seq(0,1,0.2)
  x1_labels = round(seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/100)[x_at * 100 + 1], 1) * stds[1] + means[1]
  x2_labels = round(seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/100)[x_at * 100 + 1], 1) * stds[2] + means[2]
  axis(side=1, at=x_at, labels=round(x1_labels))
  axis(side=2, at=x_at, labels=round(x2_labels))
  
  ## Add streamlines from the edges with reduced arrowheads
  num_steps = 30          # Total steps per streamline
  step_size = 0.05        # Scaling factor for movement
  skip_steps = 6          # Number of steps before placing an arrow
  num_edge_points = 6    # Number of seed points per edge
  
  # Generate starting points along edges
  x_vals = seq(0, 1, length.out = num_edge_points)
  y_vals = seq(0, 1, length.out = num_edge_points)
  
  start_points = list()
  
  # Left edge (x = 0, varying y)
  for (y in y_vals) start_points = append(start_points, list(c(0, y)))
  # Right edge (x = 1, varying y)
  for (y in y_vals) start_points = append(start_points, list(c(1, y)))
  # Bottom edge (y = 0, varying x)
  for (x in x_vals) start_points = append(start_points, list(c(x, 0)))
  # Top edge (y = 1, varying x)
  for (x in x_vals) start_points = append(start_points, list(c(x, 1)))
  
  # Loop over each streamline starting point
  add_streamlines()
    
  ## Add points
  x1_points_x = (X_[,1] - X_min[1]) / (X_max[1] - X_min[1]) 
  x2_points_x = (X_[,2] - X_min[2]) / (X_max[2] - X_min[2]) 
  points(x1_points_x, x2_points_x, pch=X_[,5]+1, col=c('black','red')[X_[,5]+1])
  
  ## Contours
  contour((predictions_N[[k]]), levels=0, add=T, col='royalblue', lwd=2)
  contour((predictions_Z[[k]]), levels=0, add=T, col='red', lwd=2)  
    
}
#
par(mfrow=c(1,1))

#
###

###################################
## COMBINED ANALYSIS - VERSION 2 ##
###################################

## Graphical parameters
breaks = c(1,2,3,4,5)

## Figure
par(mfrow=c(3,3), bty='o')
#
k = 1
## Plot
img = (predictions_N[[k]] > 0) + (predictions_Z[[k]] < 0) * 2
image(img*0, col='white', xaxt='n', yaxt='n', xlab='Density (fish/tank)', ylab='Body Mass (mg)')
# image(img*0, col=hcl.colors(length(breaks)-1, "RdYlBu", rev=FALSE, alpha = 0.3), xaxt='n', yaxt='n', xlab='Log Density', ylab='Log Body Mass')
# image(img, col=hcl.colors(length(breaks)-1, "RdBu", rev=FALSE, alpha = 0.5), xaxt='n', yaxt='n', xlab='Log Density', ylab='Log Body Mass')
# contour((predictions_N[[k]]) * (predictions_Z[[k]]), levels=0, add=T, col='red', lwd=4)  

# ## Add points  
# pch_vector = c()
# for (x in 1:nrow(img)){
#   for(y in 1:ncol(img)){
#     points(x/nrow(img), y/ncol(img), pch=img[x,y])
#   }
# }

## Axes
x_at = seq(0,1,0.2)
x1_labels = round(seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/100)[x_at * 100 + 1], 1) * stds[1] + means[1]
x2_labels = round(seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/100)[x_at * 100 + 1], 1) * stds[2] + means[2]
axis(side=1, at=x_at, labels=round(x1_labels))
axis(side=2, at=x_at, labels=round(x2_labels))

## Add points
x1_points_x = (X_[,1] - X_min[1]) / (X_max[1] - X_min[1]) 
x2_points_x = (X_[,2] - X_min[2]) / (X_max[2] - X_min[2]) 
points(x1_points_x, x2_points_x, pch=X_[,5]+1, col=c('black','red')[X_[,5]+1])

## Contours
i = 1 # two levels of fishing (0,1)
j = 2 # three levels of winter temperature (low, mid, high)
k = 1:3 # three levels of summer temperature (low, mid, high)
m = 1
for (l in c((i-1)*9 + (j-1)*3 + k*1)){
  print(l)
  contour((predictions_N[[l]]), levels=0, add=T, col='royalblue', lwd=2, lty=c(3,2,1)[m])
  contour((predictions_Z[[l]]), levels=0, add=T, col='red', lwd=2, lty=c(3,2,1)[m])  
  m = m + 1
}

## Plot 2
k = 1
img = (predictions_N[[k]] > 0) + (predictions_Z[[k]] < 0) * 2
image(img*0, col='white', xaxt='n', yaxt='n', xlab='Density (fish/tank)', ylab='Body Mass (mg)')

## Axes
x_at = seq(0,1,0.2)
x1_labels = round(seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/100)[x_at * 100 + 1], 1) * stds[1] + means[1]
x2_labels = round(seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/100)[x_at * 100 + 1], 1) * stds[2] + means[2]
axis(side=1, at=x_at, labels=round(x1_labels))
axis(side=2, at=x_at, labels=round(x2_labels))

## Add points
x1_points_x = (X_[,1] - X_min[1]) / (X_max[1] - X_min[1]) 
x2_points_x = (X_[,2] - X_min[2]) / (X_max[2] - X_min[2]) 
points(x1_points_x, x2_points_x, pch=X_[,5]+1, col=c('black','red')[X_[,5]+1])

## Contours
i = 1 # two levels of fishing (0,1)
j = 1:3 # three levels of winter temperature (low, mid, high)
k = 2 # three levels of summer temperature (low, mid, high)
m = 1
for (l in c((i-1)*9 + (j-1)*3 + k*1)){
  print(l)
  contour((predictions_N[[l]]), levels=0, add=T, col='royalblue', lwd=2, lty=c(3,2,1)[m])
  contour((predictions_Z[[l]]), levels=0, add=T, col='red', lwd=2, lty=c(3,2,1)[m])  
  m = m + 1
}

## Plot 3
k = 1
img = (predictions_N[[k]] > 0) + (predictions_Z[[k]] < 0) * 2
image(img*0, col='white', xaxt='n', yaxt='n', xlab='Density (fish/tank)', ylab='Body Mass (mg)')

## Axes
x_at = seq(0,1,0.2)
x1_labels = round(seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/100)[x_at * 100 + 1], 1) * stds[1] + means[1]
x2_labels = round(seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/100)[x_at * 100 + 1], 1) * stds[2] + means[2]
axis(side=1, at=x_at, labels=round(x1_labels))
axis(side=2, at=x_at, labels=round(x2_labels))

## Add points
x1_points_x = (X_[,1] - X_min[1]) / (X_max[1] - X_min[1]) 
x2_points_x = (X_[,2] - X_min[2]) / (X_max[2] - X_min[2]) 
points(x1_points_x, x2_points_x, pch=X_[,5]+1, col=c('black','red')[X_[,5]+1])

## Contours
i = 1:2 # two levels of fishing (0,1)
j = 2 # three levels of winter temperature (low, mid, high)
k = 2 # three levels of summer temperature (low, mid, high)
m = 1
for (l in c((i-1)*9 + (j-1)*3 + k*1)){
  print(l)
  contour((predictions_N[[l]]), levels=0, add=T, col='royalblue', lwd=2, lty=c(3,1)[m])
  contour((predictions_Z[[l]]), levels=0, add=T, col='red', lwd=2, lty=c(3,1)[m])  
  m = m + 1
}

# ## Add streamlines from the edges with reduced arrowheads
# num_steps = 30          # Total steps per streamline
# step_size = 0.05        # Scaling factor for movement
# skip_steps = 6          # Number of steps before placing an arrow
# num_edge_points = 6    # Number of seed points per edge
# 
# # Generate starting points along edges
# x_vals = seq(0, 1, length.out = num_edge_points)
# y_vals = seq(0, 1, length.out = num_edge_points)
# 
# start_points = list()
# 
# # Left edge (x = 0, varying y)
# for (y in y_vals) start_points = append(start_points, list(c(0, y)))
# # Right edge (x = 1, varying y)
# for (y in y_vals) start_points = append(start_points, list(c(1, y)))
# # Bottom edge (y = 0, varying x)
# for (x in x_vals) start_points = append(start_points, list(c(x, 0)))
# # Top edge (y = 1, varying x)
# for (x in x_vals) start_points = append(start_points, list(c(x, 1)))
# 
# Loop over each streamline starting point
# add_streamlines()

# index = NULL 
# for (i in 1:2) for (j in 1:3) for (k in 1:3) index = rbind(index, c(i,j,k))

#
par(mfrow=c(1,1))

#
###

###################################
## COMBINED ANALYSIS - VERSION 3 ##
###################################

## Graphical parameters
breaks = c(1,2,3,4,5)

index = NULL
i = 1; j = 1:3; k = 2
index = c(index, (i-1)*9 + (j-1)*3 + k)
i = 2; j = 1:3; k = 2
index = c(index, (i-1)*9 + (j-1)*3 + k)

## Figure
par(mfrow=c(3,3), bty='o')
#
for(k in index){
  
  ## Plot
  img = (predictions_N[[k]] > 0) + (predictions_Z[[k]] < 0) * 2
  image(img, col=hcl.colors(length(breaks)-1, "RdYlBu", rev=FALSE, alpha = 0.3), xaxt='n', yaxt='n', xlab='Density (fish/tank)', ylab='Body Mass (mg)')
  # image(img, col=hcl.colors(length(breaks)-1, "RdBu", rev=FALSE, alpha = 0.5), xaxt='n', yaxt='n', xlab='Log Density', ylab='Log Body Mass')
  # contour((predictions_N[[k]]) * (predictions_Z[[k]]), levels=0, add=T, col='red', lwd=4)  
  
  ## Axes
  x_at = seq(0,1,0.2)
  x1_labels = round(seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/100)[x_at * 100 + 1], 1) * stds[1] + means[1]
  x2_labels = round(seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/100)[x_at * 100 + 1], 1) * stds[2] + means[2]
  axis(side=1, at=x_at, labels=round(x1_labels))
  axis(side=2, at=x_at, labels=round(x2_labels))
  
  ## Add streamlines from the edges with reduced arrowheads
  num_steps = 30          # Total steps per streamline
  step_size = 0.05        # Scaling factor for movement
  skip_steps = 6          # Number of steps before placing an arrow
  num_edge_points = 6    # Number of seed points per edge
  
  # Generate starting points along edges
  x_vals = seq(0, 1, length.out = num_edge_points)
  y_vals = seq(0, 1, length.out = num_edge_points)
  
  start_points = list()
  
  # Left edge (x = 0, varying y)
  for (y in y_vals) start_points = append(start_points, list(c(0, y)))
  # Right edge (x = 1, varying y)
  for (y in y_vals) start_points = append(start_points, list(c(1, y)))
  # Bottom edge (y = 0, varying x)
  for (x in x_vals) start_points = append(start_points, list(c(x, 0)))
  # Top edge (y = 1, varying x)
  for (x in x_vals) start_points = append(start_points, list(c(x, 1)))
  
  # Loop over each streamline starting point
  add_streamlines()
  
  ## Add points
  x1_points_x = (X_[,1] - X_min[1]) / (X_max[1] - X_min[1]) 
  x2_points_x = (X_[,2] - X_min[2]) / (X_max[2] - X_min[2]) 
  points(x1_points_x, x2_points_x, pch=X_[,5]+1, col=c('black','red')[X_[,5]+1])
  
  ## Contours
  contour((predictions_N[[k]]), levels=0, add=T, col='royalblue', lwd=2)
  contour((predictions_Z[[k]]), levels=0, add=T, col='red', lwd=2)  
  
}
#
par(mfrow=c(1,1))

#
###

dev.off()
