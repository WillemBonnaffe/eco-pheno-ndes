#####
## ##
#####

## Goal: Fit process model to time series.
## Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

## Note: Next try with simple time-differences.

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
  MTS_o_ = MTS[MTS[,1] == idx,][,-1]
  # list_o = fit_model_o(TS, sd1_o=.1, sd2_o=.1, K_o=30)
  
  ## Compute difference 
  ddt.MTS_o_ = apply(MTS_o_, 2, diff)
  MTS_o_ = MTS_o_[1:nrow(MTS_o_)-1,] # Remove last time step as used to compute difference
  ddt.MTS_o_[,1] = MTS_o_[,1]
  head(MTS_o_)
  head(ddt.MTS_o_)
  
  ## Collect objects
  MTS_o = rbind(MTS_o, cbind(idx, MTS_o_))
  ddt.MTS_o = rbind(ddt.MTS_o, cbind(idx, ddt.MTS_o_))
  
}

# ## Save results
# system("mkdir out_o/")
# write.csv(x = MTS_o, file = "out_o/MTS_o.csv")
# write.csv(x = ddt.MTS_o, file = "out_o/ddt.MTS_o.csv")

#
###

###########################
## VISUALISE TIME SERIES ##
###########################

## Dark mode
# par(bg="black", col.axis="white", col.lab="white", col.main="white", col="white")

## Plot MTS_o
num_series = length(unique(MTS_o[,1]))
num_variables = ncol(MTS_o)-2
par(mfrow=c(2,1), mar=c(3,3,1,1))
# par(mfrow=c(num_series,2), mar=c(3,3,1,1))
#
k = 2
for (k in 1:num_series)
{
  ## Time series index 
  idx = unique(MTS_o[,1])[k]
  
  ## Select time series
  s = (MTS_o[,1] == idx)
  TS_o = MTS_o[s,][,-1]
  ddt.TS_o = ddt.MTS_o[s,][,-1]
  
  ## Plot states
  colvect = rainbow(ncol(TS_o))
  for (i in 2:ncol(TS_o))
  {
    if (i == 2) plot(TS_o[,1], TS_o[,i], cex=0, ylim=c(-1,1)*3)
    lines(c(min(TS_o[,1]), max(TS_o[,1])), c(0,0), lty=2)
    x = TS_o[,1]
    y = TS_o[,i]
    # y_lo = TS_o_lo[,i]
    # y_hi = TS_o_hi[,i]
    # polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[i])
    points(TS_o[,1], TS_o[,i], col=colvect[i], pch=16)
  }
  #
  ## Plot dynamics
  colvect = rainbow(ncol(TS_o))
  for (i in 2:ncol(TS_o))
  {
    if (i == 2) plot(TS_o[,1], TS_o[,i], cex=0, ylim=c(-1,1)*3)
    lines(c(min(TS_o[,1]), max(TS_o[,1])), c(0,0), lty=2)
    x = ddt.TS_o[,1]
    y = ddt.TS_o[,i]
    # y_lo = ddt.TS_o_lo[,i]
    # y_hi = ddt.TS_o_hi[,i]
    # polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor("grey",0.2))
    lines(x, y, col=colvect[i])
  }
}
#
par(mfrow=c(1,1))

#
###

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

##################
## FORMAT MTS_o ##
##################

## format MTS_o
MTS_o = data.frame(MTS_o)
ddt.MTS_o = data.frame(ddt.MTS_o)
for(i in 3:5) MTS_o[,i] = as.numeric(MTS_o[,i])
for(i in 3:5) ddt.MTS_o[,i] = as.numeric(ddt.MTS_o[,i])
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
N = 4 # 4 explanatory variables
W_p = 10
N_p = 2 * W_p * (2 + N)
sd1_p = .1
sd2_p = .3 # N_DJK = .25; Z_DJK = .45; N_BEH = 1.0; Z_BEH = 1.25
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
  plot(x, rep(0,length(x)), ylim=c(-1,1)*4, lty=2, type="l", xlab="Time", ylab="Dynamics", main=paste("Time Series ", selected_TS_, sep=""))
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
  colvect = c("red", "green", "blue", 'magenta')
  
  ## Dynamics
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*4, lty=2, type="l", xlab="Time", ylab="Effects")
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
  lines(c(t_lb, t_lb), c(-1,1)*4, lty=2)
  lines(c(t_rb, t_rb), c(-1,1)*4, lty=2)
  #
  ## Legend
  legend("bottom", legend = c("N", "Z", "E", "F"), col=colvect, horiz=T, lty=1)
  
  ## End graph
  par(mfrow=c(1,1))
  
  #
  ###
}

################################
## COMPUTE EXPECTATION AND CI ##
################################

## Graphical parameters
par(bg="white", col.axis="black", col.lab="black", col.main="black", col="black")
# par(bg="black", col.axis="white", col.lab="white", col.main="white", col="white")
par(mfrow=c(3,3), mar=c(4.5,4.5,1,1), cex.lab=1.5)
colvect = c("red", "red", "red")
varnames = c("Log Density (S.U.)", "Log Body Length (S.U.)", "Temperature (S.U.)")
pchvect = c(1,3)

## Plot
for (F_ in c(0.0, 1.0))
{
  for (E_ in c(-1.0, 0.0, 1.0))
  {
    for (s in 1:3)
    {
      ## Get x 
      x = seq(0, 1, 0.01) * (max(X[,s]) - min(X[,s])) + min(X[,s])
      X__ = matrix(0, ncol=N, nrow=length(x))
      X__[,3] = E_
      X__[,4] = F_
      X__[,s] = x
      
      ## Expectation and credible interval
      ddt.N_p = ECI(X__, chain, ddt.n)
      ddt.N_p_interq = ECI(X__, chain, ddt.n, lb=0.25, rb=0.75)
      
      ## Dynamics
      # x = X__[,s]
      plot(x, rep(0,length(x)), ylim=c(-1,1)*2.5, lty=2, type="l", xlab=varnames[s], ylab="P.c. Population Growth")
      #
      ## Credible interval
      y = ddt.N_p$mean
      y_lo = ddt.N_p$lo
      y_hi = ddt.N_p$hi
      polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=grey(0.95))
      # polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[s],0.25))# grey(0.95))
      lines(x, y, col=colvect[s])
      #
      ## Inter-quartile
      y = ddt.N_p_interq$mean
      y_lo = ddt.N_p_interq$lo
      y_hi = ddt.N_p_interq$hi
      polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=grey(0.85))
      # polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[s],0.25))# grey(0.85))
      lines(x, y, col=colvect[s])
      #
      ## Data
      for (selected_TS_ in selected_TS)
      {
        s_ts = which((MTS_o[,1] == selected_TS_))
        X__real = X_[s_ts, ]
        Y__real = Y_[s_ts, ]
        points(X__real[,s], Y__real[,rv_idx], pch=pchvect[X__real[,4]+1], col=c('black','red')[X__real[,4]+1])
      }
      #
      ## Legend
      lines(x, rep(0,length(x)), lty=2)
      # legend("bottom", legend = varnames, col=colvect, horiz=T, lty=1, bty="n")
      legend("top", legend = paste("Temperature = ", E_, ", Fishing = ", F_), bty="n")
      legend("bottom", legend = c("Not Fished", "Fished"), horiz=T, pch=pchvect, bty="n", col=c('black','red'))
    }
  }
}
#
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
X_min = apply(X_, 2, min)
X_max = apply(X_, 2, max)
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
par(mfrow=c(2,2))
#
predictions_N = list()
k = 1
for (X_4 in c(0,1)){
  for (X_3 in c(-1.5,0.4)){
    
    ## Form mesh
    mesh = cbind(as.vector(X1), as.vector(X2), X_3, X_4)
    
    ## Predictions
    predictions = ECI(mesh, chain, ddt.n)
    predictions = predictions$mean
    predictions = matrix(predictions, nrow=res, ncol=res)
    predictions_N[[k]] = predictions
    
    ## Visualise mesh
    breaks = seq(-5,5,1.0)
    image(predictions, xaxt='n', yaxt='n', xlab='Log Density', ylab='Log Body Size', col=hcl.colors(length(breaks)-1, "PiYG", rev=FALSE, alpha = 0.75), breaks=breaks)
    contour(predictions, add=T)
    contour(predictions, add=T, levels = c(0), col = 'red', lwd=2)
    x_at = seq(0,1,0.2)
    x1_labels = round(seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/100)[x_at * 100 + 1], 2)
    x2_labels = round(seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/100)[x_at * 100 + 1], 2)
    axis(side=1, at=x_at, labels=x1_labels)
    axis(side=2, at=x_at, labels=x2_labels)
    
    ## Add points
    x1_points_x = (X_[,1] - X_min[1]) / (X_max[1] - X_min[1]) 
    x2_points_x = (X_[,2] - X_min[2]) / (X_max[2] - X_min[2]) 
    points(x1_points_x, x2_points_x, pch=X_[,4]+1)
    
    ## Iterator
    k = k + 1
    
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
N = 4 # 4 explanatory variables
W_p = 10
N_p = 2 * W_p * (2 + N)
sd1_p = .1
sd2_p = .3 # N_DJK = .25; Z_DJK = .45; N_BEH = 1.0; Z_BEH = 1.25
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

for (selected_TS_ in selected_TS)
{
  ################################
  ## COMPUTE EXPECTATION AND CI ##
  ################################
  
  ## Get time series set
  s_t = which((MTS_o[,1] == selected_TS_))
  X__ = X_[s_t, ]
  Y__ = Y_[s_t, ]
  
  ## Contributions
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
  plot(x, rep(0,length(x)), ylim=c(-1,1)*4, lty=2, type="l", xlab="Time", ylab="Dynamics", main=paste("Time Series ", selected_TS_, sep=""))
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
  colvect = c("red", "green", "blue", 'magenta')
  
  ## Dynamics
  x = 1:nrow(X__)
  plot(x, rep(0,length(x)), ylim=c(-1,1)*4, lty=2, type="l", xlab="Time", ylab="Effects")
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
  lines(c(t_lb, t_lb), c(-1,1)*4, lty=2)
  lines(c(t_rb, t_rb), c(-1,1)*4, lty=2)
  #
  ## Legend
  legend("bottom", legend = c("N", "Z", "E", "F"), col=colvect, horiz=T, lty=1)
  
  ## End graph
  par(mfrow=c(1,1))
  
  #
  ###
}

################################
## COMPUTE EXPECTATION AND CI ##
################################

## Graphical parameters
par(bg="white", col.axis="black", col.lab="black", col.main="black", col="black")
# par(bg="black", col.axis="white", col.lab="white", col.main="white", col="white")
par(mfrow=c(3,3), mar=c(4.5,4.5,1,1), cex.lab=1.5)
colvect = c("red", "red", "red")
varnames = c("Log Density (S.U.)", "Log Body Length (S.U.)", "Temperature (S.U.)")
pchvect = c(1,3)

## Plot
for (F_ in c(0.0, 1.0))
{
  for (E_ in c(-1.0, 0.0, 1.0))
  {
    for (s in 1:3)
    {
      ## Get x 
      x = seq(0, 1, 0.01) * (max(X[,s]) - min(X[,s])) + min(X[,s])
      X__ = matrix(0, ncol=N, nrow=length(x))
      X__[,3] = E_
      X__[,4] = F_
      X__[,s] = x
      
      ## Expectation and credible interval
      ddt.N_p = ECI(X__, chain, ddt.n)
      ddt.N_p_interq = ECI(X__, chain, ddt.n, lb=0.25, rb=0.75)
      
      ## Dynamics
      # x = X__[,s]
      plot(x, rep(0,length(x)), ylim=c(-1,1)*2.5, lty=2, type="l", xlab=varnames[s], ylab="Phenotypic Change")
      #
      ## Credible interval
      y = ddt.N_p$mean
      y_lo = ddt.N_p$lo
      y_hi = ddt.N_p$hi
      polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=grey(0.95))
      # polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[s],0.25))# grey(0.95))
      lines(x, y, col=colvect[s])
      #
      ## Inter-quartile
      y = ddt.N_p_interq$mean
      y_lo = ddt.N_p_interq$lo
      y_hi = ddt.N_p_interq$hi
      polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=grey(0.85))
      # polygon(x=c(x,rev(x)), y=c(y_lo, rev(y_hi)), border=NA, col=adjustcolor(colvect[s],0.25))# grey(0.85))
      lines(x, y, col=colvect[s])
      #
      ## Data
      for (selected_TS_ in selected_TS)
      {
        s_ts = which((MTS_o[,1] == selected_TS_))
        X__real = X_[s_ts, ]
        Y__real = Y_[s_ts, ]
        points(X__real[,s], Y__real[,rv_idx], pch=pchvect[X__real[,4]+1], col=c('black','red')[X__real[,4]+1])
      }
      #
      ## Legend
      lines(x, rep(0,length(x)), lty=2)
      # legend("bottom", legend = varnames, col=colvect, horiz=T, lty=1, bty="n")
      legend("top", legend = paste("Temperature = ", E_, ", Fishing = ", F_), bty="n")
      legend("bottom", legend = c("Not Fished", "Fished"), horiz=T, pch=pchvect, bty="n", col=c('black','red'))
    }
  }
}
#
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
X_min = apply(X_, 2, min)
X_max = apply(X_, 2, max)
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
par(mfrow=c(2,2))
#
predictions_Z = list()
k = 1
for (X_4 in c(0,1)){
  for (X_3 in c(-1.5,0.4)){
    
    ## Form mesh
    mesh = cbind(as.vector(X1), as.vector(X2), X_3, X_4)
    
    ## Predictions
    predictions = ECI(mesh, chain, ddt.n)
    predictions = predictions$mean
    predictions = matrix(predictions, nrow=res, ncol=res)
    predictions_Z[[k]] = predictions
    
    ## Visualise mesh
    breaks = seq(-5,5,1.0)
    image(predictions, xaxt='n', yaxt='n', xlab='Log Density', ylab='Log Body Size', col=hcl.colors(length(breaks)-1, "PiYG", rev=FALSE, alpha = 0.75), breaks=breaks)
    contour(predictions, add=T)
    contour(predictions, add=T, levels = c(0), col = 'red', lwd=2)
    
    ## Axes
    x_at = seq(0,1,0.2)
    x1_labels = round(seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/100)[x_at * 100 + 1], 2)
    x2_labels = round(seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/100)[x_at * 100 + 1], 2)
    axis(side=1, at=x_at, labels=x1_labels)
    axis(side=2, at=x_at, labels=x2_labels)
    
    ## Add points
    x1_points_x = (X_[,1] - X_min[1]) / (X_max[1] - X_min[1]) 
    x2_points_x = (X_[,2] - X_min[2]) / (X_max[2] - X_min[2]) 
    points(x1_points_x, x2_points_x, pch=X_[,4]+1)
    
    ## Iterator
    k = k + 1
    
  }
}
#
par(mfrow=c(1,1)) 

#
###

#######################
## COMBINED ANALYSIS ##
#######################

## Graphical parameters
breaks = c(1,2,3,4,5)

## Figure
par(mfrow=c(2,2))
#
for(k in 1:4){
  
  ## Plot
  img = (predictions_N[[k]] < 0) + (predictions_Z[[k]] > 0) * 2
  image(img, col=hcl.colors(length(breaks)-1, "PiYG", rev=FALSE, alpha = 0.5), xaxt='n', yaxt='n', xlab='Log Density', ylab='Log Body Size')
  # image(img, col=hcl.colors(length(breaks)-1, "RdBu", rev=FALSE, alpha = 0.5), xaxt='n', yaxt='n', xlab='Log Density', ylab='Log Body Size')
  # contour((predictions_N[[k]]) * (predictions_Z[[k]]), levels=0, add=T, col='red', lwd=4)  
  
  ## Contours
  contour((predictions_N[[k]]), levels=0, add=T, col='green', lwd=4)  
  contour((predictions_Z[[k]]), levels=0, add=T, col='red', lwd=4)  
  
  ## Axes
  x_at = seq(0,1,0.2)
  x1_labels = round(seq(X_min[1], X_max[1], (X_max[1] - X_min[1])/100)[x_at * 100 + 1], 1)
  x2_labels = round(seq(X_min[2], X_max[2], (X_max[2] - X_min[2])/100)[x_at * 100 + 1], 1)
  axis(side=1, at=x_at, labels=x1_labels)
  axis(side=2, at=x_at, labels=x2_labels)
  
  ## Add points
  x1_points_x = (X_[,1] - X_min[1]) / (X_max[1] - X_min[1]) 
  x2_points_x = (X_[,2] - X_min[2]) / (X_max[2] - X_min[2]) 
  points(x1_points_x, x2_points_x, pch=X_[,4]+1)
  
  ## Add streamlines from the edges with reduced arrowheads
  num_steps = 20          # Total steps per streamline
  step_size = 0.05        # Scaling factor for movement
  skip_steps = 3          # Number of steps before placing an arrow
  num_edge_points = 15    # Number of seed points per edge
  
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
  for (point in start_points) {
    
    # Initialize position
    x_pos = point[1]
    y_pos = point[2]
    
    for (step in 1:num_steps) {
      # Find the closest grid index
      i = round(x_pos * (res - 1)) + 1
      j = round(y_pos * (res - 1)) + 1
      
      # Ensure indices are within valid range
      if (i < 1 || i > res || j < 1 || j > res) break
      
      # Get vector field values
      dx = predictions_N[[k]][i, j] * step_size
      dy = predictions_Z[[k]][i, j] * step_size
      
      # Compute new position
      x_new = x_pos + dx
      y_new = y_pos + dy
      
      # Ensure new position stays within bounds
      if (x_new < 0 || x_new > 1 || y_new < 0 || y_new > 1) break
      
      # Draw an arrow every `skip_steps` steps, otherwise just a line
      if (step %% skip_steps == 0) {
        arrows(x_pos, y_pos, x_new, y_new, length=0.05, col=adjustcolor("gray20", alpha=1.0), lwd=1.0)
      } else {
        segments(x_pos, y_pos, x_new, y_new, col=adjustcolor("gray20", alpha=1.0), lwd=1.0)
      }
      
      # Update position
      x_pos = x_new
      y_pos = y_new
      
    }
    
    ## Contours
    contour((predictions_N[[k]]), levels=0, add=T, col='green', lwd=4)  
    contour((predictions_Z[[k]]), levels=0, add=T, col='red', lwd=4)  
    
  }
}
#
par(mfrow=c(1,1))

#
###