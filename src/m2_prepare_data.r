#####
## ##
#####

## Goal: 
## Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

##############
## INITIATE ##
##############

## Local imports
source("f_betteRplots.r")
source("f_utils.r")

#
###

##############
## INITIATE ##
##############

## Goal: load data, functions

## User defined parameters
timeSeriesIds = c("A","B","C","D","E","F","G","H","I","J","K","L")
selected_columns = c("t","N","Z_mean","temp_mean_summer","temp_mean_winter")
column_names = c("t","N","Z","E1","E2")

## Compute properties of time series
MTS = NULL
# n_time_steps = NULL
for (timeSeriesId in timeSeriesIds)
{
  ## Load data
  TS = read.table(paste("data/TS_", timeSeriesId, ".csv", sep=""),sep=",",header=T)[selected_columns]
  
  ## Subset spring sample only
  s = seq(1+1, nrow(TS)+1, 2)
  TS = TS[s,]
  
  ## Remove first time step (i.e. introduction)
  TS = TS[TS$t!=0.5,]

  ## Add time series Id
  MTS_ = cbind(timeSeriesId, TS)
  
  ## Collect
  MTS = rbind(MTS, MTS_)
}

# ## Set minimum value
# th = .1
# s = c(3,4)
# MTS[,s] = MTS[,s] * (MTS[,s] > th) + (MTS[,s] < th) * th
# 
# ## Transform
# s = c(3,4)
# MTS[,s] = log(MTS[,s])
# 
# ## Compute mean
# s = -c(1,2)
# means = apply(MTS[,s], 2, mean, na.rm=T)
# stds = apply(MTS[,s], 2, sd, na.rm=T)
# 
# ## Standardise
# MTS_ = MTS
# s = -c(1,2)
# MTS[,s] = t((t(MTS[,s])-means)/stds)
# plot(MTS[,5], MTS_[,5]) # Check

# # Retain only odd time steps
# s = NULL
# selected = c("0.5","1.5","2.5","3.5","4.5","5.5")
# for (selected_ in selected)
# {
#   s_ = which((MTS[,2] == selected_))
#   s = c(s, s_)
# }
# MTS = MTS[s,]
# MTS = MTS[order(MTS[,1]),]

## Format
head(MTS)

## Store
write.csv(x = MTS, file = "data/MTS_all.csv", row.names = F)

#
###

###########################
## VISUALISE TIME SERIES ##
###########################

## Parameters
num_col = ncol(MTS)

## Dark mode
# par(bg="black", col.axis="white", col.lab="white", col.main="white", col="white")
par(bg="white", col.axis="black", col.lab="black", col.main="black", col="black")

## Plot time series (V1)
par(mfrow=c(num_col-1,num_col-1))
colvect = rainbow(length(unique(MTS[,1])))
for (i in 2:num_col)
{
  for (j in 2:num_col)
  {
    if (i == j)
    {
      plot(-1:1,-1:1, cex=0, xaxt="n", yaxt="n", xlab='', ylab='')
      text(0,0,paste(colnames(MTS)[i]), cex=2)
    } else
    {
      plot(MTS[,j], MTS[,i], cex=0, xlab='', ylab='')
      for (k in 1:length(unique(MTS[,1])))
      {
        s = MTS[,1] == unique(MTS[,1])[k]
        lines(MTS[s,j], MTS[s,i], col=colvect[k], type="b", pch=16)
      }
    }
  }
}
par(mfrow=c(1,1))

## Plot time series (V1)
num_variables = ncol(MTS) - 2
num_series = length(unique(MTS[,1]))
colvect = rainbow(length(unique(MTS[,1])))
varnames = c('Log Density', 'Log Body Length', 'Summer Mean Temperature', 'Winter Mean Temperature')
par(mfrow=c(num_variables,1))
for (i in 1:num_variables)
{
  plot(MTS[,2], MTS[,i+2], cex=0, xlab="Time (Years)", ylab=varnames[i])
  for (j in 1:num_series)
  {
    s = MTS[,1] == unique(MTS[,1])[j]
    lines(MTS[s,2], MTS[s,i+2], col=colvect[j], type="b", pch=16)
  }
  legend('top', horiz=T, legend=timeSeriesIds, lty=1, col=colvect, bty="n")
}
par(mfrow=c(1,1))

#
###