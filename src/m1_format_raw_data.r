#####
## ##
#####

## Goal:
## Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

##############
## INITIATE ##
##############

dat = read.table("raw/Masses_all.csv", sep=",", header=T)
head(dat)
nrow(dat)
ncol(dat)

#
###

#################
## FORMAT DATA ##
#################

## get sampling event index by year
yearly_sampling_event = rep(0,nrow(dat))
for (i in 2:4) yearly_sampling_event[which(dat$Month==i)] = 1
for (i in 10:11) yearly_sampling_event[which(dat$Month==i)] = 2
plot(yearly_sampling_event)
sampling_event_by_month = dat$Year*12 + yearly_sampling_event*6 - min(dat$Year)*12
sampling_event_by_year = sampling_event_by_month/12
plot(sampling_event_by_year)
dat$sampling_event_by_year = sampling_event_by_year

## initate aggregated data
aggregated_data = NULL

## for each tank
unique_tanks = unique(dat$Tank)
for (unique_tank in unique_tanks)
{
  ## check
  print(unique_tank)
  dat_ = dat[dat$Tank == unique_tank,]
  
  ## initiate level _
  aggregated_data_ = NULL
  
  ## for each sampling event
  unique_sampling_events = unique(dat_$sampling_event_by_year)
  for (unique_sampling_event in unique_sampling_events)
  {
    ## check
    print(unique_sampling_event)
    dat__ = dat_[dat_$sampling_event_by_year==unique_sampling_event,]
    n__ = nrow(dat__)
    mean_mass__ = mean(dat__$Mass, na.rm=T)
    sd_mass__ = sd(dat__$Mass, na.rm=T)
    
    ## initiate level __
    aggregated_data__ = cbind(unique_tank, unique_sampling_event, n__, mean_mass__, sd_mass__)
    
    ## aggregate level __
    aggregated_data_ = rbind(aggregated_data_, aggregated_data__)
  }
  ## aggregate level _
  aggregated_data = rbind(aggregated_data, aggregated_data_)
}

## format 
aggregated_data = data.frame(aggregated_data)
colnames(aggregated_data) = c("tank", "t", "N", "Z_mean", "Z_sd")

#
###

####################
## VISUALISATIONS ##
####################

##
TS = aggregated_data

## visualisations of the time series
s = TS$tank == "A"
plot(TS[s,-1])

## visualise the time series in parallel
par(mfrow=c(3,3))
for (tank in unique(TS$tank))
{
  s = TS$tank == tank
  plot(TS[s,]$t, TS[s,]$N, type="l", xlab="Time", ylab="N", main=paste(tank))
  plot(TS[s,]$t, TS[s,]$Z_mean, type="l", xlab="Time", ylab="Z Mean")
  plot(TS[s,]$t, TS[s,]$Z_sd, type="l", xlab="Time", ylab="Z St. Deviation")
}
par(mfrow=c(1,1))

## visualise the time series stacked
par(mfrow=c(3,3))
#
## population size
k = 2
plot(TS$t, TS$N, xlab="Time", ylab="N", cex=0)
for (tank in unique(TS$tank))
{
  s = TS$tank == tank
  lines(TS[s,]$t, TS[s,]$N, type="l", col=k)
  k = k + 1
}
#
## mean phenotype
k = 2
plot(TS$t, TS$Z_mean, xlab="Time", ylab="Z Mean", cex=0)
for (tank in unique(TS$tank))
{
  s = TS$tank == tank
  lines(TS[s,]$t, TS[s,]$Z_mean, type="l", col=k)
  k = k + 1
}
#
## phenotype standard deviation
k = 2
plot(TS$t, TS$Z_sd, xlab="Time", ylab="Z St. Deviation", cex=0)
for (tank in unique(TS$tank))
{
  s = TS$tank == tank
  lines(TS[s,]$t, TS[s,]$Z_sd, type="l", col=k)
  k = k + 1
}
par(mfrow=c(1,1))

#
###

# ##########
# ## SAVE ##
# ##########
# 
# ## save time series of each tank
# for (tank in unique_tanks)
# {
#   ## for each tank
#   print(tank)
#   s = TS$tank == tank
#   TS_ = TS[s,-1]
#   
#   ## format save
#   TS_ = cbind(rep(tank, nrow(TS_)), TS_)
#   colnames(TS_) = c("tank", "t", "N", "Z_mean", "Z_sd")
#   write.table(x = TS_, file = paste("TS_", tank, ".csv",sep=""), sep = ",", row.names = FALSE)
# }
# 
# #
# ###

#################################
## ADD ENVIRONMENTAL VARIABLES ##
#################################

## load data
env_dat = read.table("raw/Meteo_France_climate_at_enclosures.txt", sep="\t", header=T)
# sel_dat = env_dat[,7:18] # selection differential
head(env_dat)
env_dat = env_dat[,c("Year", "Mean.summer.temp", "Mean.winter.temp")] # only keep environmental variables
nrow(env_dat)
plot(env_dat$Year, env_dat$Mean.summer.temp, type='l', col="salmon", ylim=c(0,20))
lines(env_dat$Year, env_dat$Mean.winter.temp, type='l', col="blue")

## Add missing rows (summer census)
env_dat_ = NULL
# env_dat_ = rbind(env_dat_, as.numeric(rep(NA, ncol(env_dat))))
# env_dat_ = rbind(env_dat_, as.numeric(rep(NA, ncol(env_dat))))
for (i in 1:nrow(env_dat))
{
  env_dat_ = rbind(env_dat_, as.numeric(rep(NA, ncol(env_dat))))
  env_dat_ = rbind(env_dat_, as.numeric(env_dat[i,]))
}
env_dat_ = rbind(env_dat_, as.numeric(rep(NA, ncol(env_dat))))
env_dat_ = rbind(env_dat_, as.numeric(rep(NA, ncol(env_dat))))
env_dat_ = data.frame(env_dat_)
colnames(env_dat_) = colnames(env_dat)
env_dat = env_dat_
nrow(env_dat)

## save time series of each tank
for (tank in unique_tanks)
{
  ## for each tank
  print(tank)
  s = TS$tank == tank
  TS_ = TS[s,-1]
  
  ## concatenate environmental variables
  TS_ = cbind(TS_, env_dat)
  
  ## format and save
  TS_ = cbind(rep(tank, nrow(TS_)), TS_)
  colnames(TS_) = c("tank", "t", "N", "Z_mean", "Z_sd", "year", "temp_mean_summer",  "temp_mean_winter")
  write.table(x = TS_, file = paste("data/TS_", tank, ".csv",sep=""), sep = ",", row.names = FALSE)
}

#
###