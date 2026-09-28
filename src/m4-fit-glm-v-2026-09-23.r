#####
## ##
#####

## Goal: Fit process model to time series.
## Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

## Update log:
## 2025-09-02 - Added collection of table of effects for ∆N(t) and ∆Z(t).
## 2026-02-21 - Added de-transformation of variables for visualisations.
## 2026-09-23 - Compare results produced by simple multiple linear regression

pdf(file = "figures-2026-09-23.pdf", width = 12, height = 12)

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

####################################
## PREPARE DATA FOR LINEAR MODELS ##
####################################

dat <- data.frame(X_)
dat$ddt.pop <- Y_$N
dat$ddt.phe <- Y_$Z_mean
colnames(dat)

#
###

####################
## LINEAR MODEL 1 ##
####################

## Backward step-wise model simplification
LM1.0 <- lm(ddt.pop ~ N + Z_mean + temp_mean_summer * fishing_status + temp_mean_winter * fishing_status, data = dat)
summary(LM1.0)
LM1.1 <- lm(ddt.pop ~ N + Z_mean + temp_mean_summer * fishing_status + temp_mean_winter, data = dat)
summary(LM1.1)
LM1.2 <- lm(ddt.pop ~ N + temp_mean_summer * fishing_status + temp_mean_winter, data = dat)
summary(LM1.2)

## Check AIC
AIC(LM1.0, k=2)
AIC(LM1.1, k=2)
AIC(LM1.2, k=2)

## Check residuals
plot(LM1.2)

## Simplification from full
LM1.full.0 <- lm(ddt.pop ~ N * fishing_status + Z_mean * fishing_status + temp_mean_summer * fishing_status + temp_mean_winter * fishing_status, data = dat)
summary(LM1.full.0)
drop1(LM1.full.0)
LM1.full.1 <- lm(ddt.pop ~ N * fishing_status + Z_mean * fishing_status + temp_mean_summer + temp_mean_winter * fishing_status, data = dat)
summary(LM1.full.1)
drop1(LM1.full.1)
LM1.full.2 <- lm(ddt.pop ~ N * fishing_status + Z_mean + temp_mean_summer + temp_mean_winter * fishing_status, data = dat)
summary(LM1.full.2)
drop1(LM1.full.2)
LM1.full.3 <- lm(ddt.pop ~ N * fishing_status + Z_mean + temp_mean_summer + temp_mean_winter, data = dat)
summary(LM1.full.3)
drop1(LM1.full.3)
LM1.full.4 <- lm(ddt.pop ~ N * fishing_status + temp_mean_summer + temp_mean_winter, data = dat)
summary(LM1.full.4)
drop1(LM1.full.4)

## Check AIC
AIC(LM1.full.0, k=2)
AIC(LM1.full.1, k=2)
AIC(LM1.full.2, k=2)
AIC(LM1.full.3, k=2)
AIC(LM1.full.4, k=2)

## Check changes in residuals
anova(LM1.full.1,LM1.full.0)
anova(LM1.full.2,LM1.full.1)
anova(LM1.full.4,LM1.full.3)

#
###

###########################
## VISUALISE PREDICTIONS ##
###########################

plot_effect <- function(var, label, model, response) {
  
  x <- dat[[var]]
  
  plot(
    x[dat$fishing_status == 0],
    response[dat$fishing_status == 0],
    pch = 1, col = "black",
    xlab = label, ylab = "Population growth",
    xlim = range(x, na.rm = TRUE),
    ylim = range(dat$ddt.pop, na.rm = TRUE),
  )
  
  points(
    x[dat$fishing_status == 1],
    response[dat$fishing_status == 1],
    pch = 2, col = "red"
  )
  
  newdat <- data.frame(
    N = mean(dat$N, na.rm = TRUE),
    Z_mean = mean(dat$Z_mean, na.rm = TRUE),
    temp_mean_summer = mean(dat$temp_mean_summer, na.rm = TRUE),
    temp_mean_winter = mean(dat$temp_mean_winter, na.rm = TRUE),
    fishing_status = rep(c(0, 1), each = 100)
  )
  
  newdat[[var]] <- rep(
    seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100),
    2
  )
  
  pred <- predict(model, newdata = newdat)
  
  lines(newdat[[var]][newdat$fishing_status == 0],
        pred[newdat$fishing_status == 0],
        col = "black", lwd = 2)
  
  lines(newdat[[var]][newdat$fishing_status == 1],
        pred[newdat$fishing_status == 1],
        col = "red", lwd = 2)
}

par(mfrow = c(2, 2))
plot_effect("N", "Density", LM1.full.4, dat$ddt.pop)
lines(c(-10,10),c(0,0),lty=2)
plot_effect("Z_mean", "Mean phenotype", LM1.full.4, dat$ddt.pop)
lines(c(-10,10),c(0,0),lty=2)
plot_effect("temp_mean_summer", "Summer temperature", LM1.full.4, dat$ddt.pop)
lines(c(-10,10),c(0,0),lty=2)
plot_effect("temp_mean_winter", "Winter temperature", LM1.full.4, dat$ddt.pop)
lines(c(-10,10),c(0,0),lty=2)

legend(
  "topright",
  legend = c("Non-harvested", "Harvested"),
  pch = c(1, 17),
  col = c("black", "red"),
  bty = "n"
)

#
###

####################
## LINEAR MODEL 2 ##
####################

## Backward step-wise model simplification
LM2.0 <- lm(ddt.phe ~ N + Z_mean + temp_mean_summer * fishing_status + temp_mean_winter * fishing_status, data = dat)
summary(LM2.0)
LM2.1 <- lm(ddt.phe ~ N + Z_mean + temp_mean_summer * fishing_status + temp_mean_winter, data = dat)
summary(LM2.1)
LM2.2 <- lm(ddt.phe ~ N + Z_mean + temp_mean_summer + temp_mean_winter + fishing_status, data = dat)
summary(LM2.2)
drop1(LM2.2)
LM2.3 <- lm(ddt.phe ~ Z_mean + temp_mean_summer + temp_mean_winter + fishing_status, data = dat)
summary(LM2.3)
LM2.4 <- lm(ddt.phe ~ Z_mean + temp_mean_winter + fishing_status, data = dat)
summary(LM2.4)
LM2.5 <- lm(ddt.phe ~ Z_mean + temp_mean_winter, data = dat)
summary(LM2.5)

## Check AIC
AIC(LM2.0, k=2)
AIC(LM2.1, k=2)
AIC(LM2.2, k=2)
AIC(LM2.3, k=2)
AIC(LM2.4, k=2)
AIC(LM2.5, k=2)

## Check residuals
plot(LM2.5)

## Simplification from full
LM2.full.0 <- lm(ddt.phe ~ N * fishing_status + Z_mean * fishing_status + temp_mean_summer * fishing_status + temp_mean_winter * fishing_status, data = dat)
summary(LM2.full.0)
drop1(LM2.full.0)
LM2.full.1 <- lm(ddt.phe ~ N + Z_mean * fishing_status + temp_mean_summer * fishing_status + temp_mean_winter * fishing_status, data = dat)
summary(LM2.full.1)
drop1(LM2.full.1)
LM2.full.2 <- lm(ddt.phe ~ N + Z_mean * fishing_status + temp_mean_summer + temp_mean_winter * fishing_status, data = dat)
summary(LM2.full.2)
drop1(LM2.full.2)
LM2.full.3 <- lm(ddt.phe ~ N + Z_mean * fishing_status + temp_mean_summer + temp_mean_winter, data = dat)
summary(LM2.full.3)
drop1(LM2.full.3)
LM2.full.4 <- lm(ddt.phe ~ Z_mean * fishing_status + temp_mean_summer + temp_mean_winter, data = dat)
summary(LM2.full.4)
drop1(LM2.full.4)
LM2.full.5 <- lm(ddt.phe ~ Z_mean * fishing_status + temp_mean_winter, data = dat)
summary(LM2.full.5)
drop1(LM2.full.5)

## Check AIC
AIC(LM2.full.0, k=2)
AIC(LM2.full.1, k=2)
AIC(LM2.full.2, k=2)
AIC(LM2.full.3, k=2)
AIC(LM2.full.4, k=2)
AIC(LM2.full.5, k=2)

## Check changes in residuals
anova(LM2.full.1,LM2.full.0)
anova(LM2.full.2,LM2.full.1)
anova(LM2.full.4,LM2.full.3)
anova(LM2.full.5,LM2.full.4)

#
###

###########################
## VISUALISE PREDICTIONS ##
###########################

par(mfrow = c(2, 2))
plot_effect("N", "Density", LM2.full.5, dat$ddt.phe)
lines(c(-10,10),c(0,0),lty=2)
plot_effect("Z_mean", "Mean phenotype", LM2.full.5, dat$ddt.phe)
lines(c(-10,10),c(0,0),lty=2)
plot_effect("temp_mean_summer", "Summer temperature", LM2.full.5, dat$ddt.phe)
lines(c(-10,10),c(0,0),lty=2)
plot_effect("temp_mean_winter", "Winter temperature", LM2.full.5, dat$ddt.phe)
lines(c(-10,10),c(0,0),lty=2)

legend(
  "topright",
  legend = c("Non-harvested", "Harvested"),
  pch = c(1, 17),
  col = c("black", "red"),
  bty = "n"
)

#
###

########################################
## CHECK TEMPORAL AUTOCORRELATION ##
########################################

## Residuals are checked separately for each time series
time_series <- split(seq_len(nrow(dat)), MTS_o$idx)

check_residual_autocorrelation <- function(model, model_name) {
  do.call(rbind, lapply(names(time_series), function(id) {
    residuals <- resid(model)[time_series[[id]]]
    acf_result <- acf(residuals, plot = FALSE)
    test_result <- Box.test(residuals, lag = 3, type = "Ljung-Box")

    data.frame(
      Model = model_name,
      Time_series = id,
      Lag_1_ACF = acf_result$acf[2],
      Ljung_Box_p_value = test_result$p.value
    )
  }))
}

residual_autocorrelation <- rbind(
  check_residual_autocorrelation(LM1.full.4, "LM1_population"),
  check_residual_autocorrelation(LM2.full.5, "LM2_phenotype")
)

write.csv(
  residual_autocorrelation,
  "table_residual_temporal_autocorrelation.csv",
  row.names = FALSE
)

## ACF plots for each time series
pdf("figures_residual_temporal_autocorrelation.pdf", width = 10, height = 8)
par(mfrow = c(3, 4))

for (id in names(time_series)) {
  acf(resid(LM1.full.4)[time_series[[id]]],
      main = paste("LM1 residuals:", id))
}

for (id in names(time_series)) {
  acf(resid(LM2.full.5)[time_series[[id]]],
      main = paste("LM2 residuals:", id))
}

dev.off()

#
###

###########################
## EXPORT MODEL TABLES ##
###########################

## Parameter estimates for the final model in each model sequence
make_parameter_table <- function(model) {
  x <- as.data.frame(coef(summary(model)))
  x$Parameter <- rownames(x)
  rownames(x) <- NULL
  x <- x[, c("Parameter", "Estimate", "Std. Error", "Pr(>|t|)")]
  colnames(x) <- c("Parameter", "Estimate", "Std_Error", "Significance")
  x
}

## Changes between successive models in a simplification sequence
make_simplification_table <- function(models, deleted_terms) {
  out <- vector("list", length(deleted_terms))

  for (i in seq_along(deleted_terms)) {
    full_model <- models[[i]]
    reduced_model <- models[[i + 1]]
    a <- anova(reduced_model, full_model)

    out[[i]] <- data.frame(
      Initial_full_model = names(models)[i],
      Deleted_term = deleted_terms[i],
      Delta_AIC = as.numeric(AIC(reduced_model)) - as.numeric(AIC(full_model)),
      Delta_Sum_of_Sq = a$`Sum of Sq`[2],
      F_statistic = a$F[2],
      Significance = a$`Pr(>F)`[2]
    )
  }

  do.call(rbind, out)
}

## Model 1: population growth
LM1_parameters <- make_parameter_table(LM1.full.4)
LM1_simplification <- make_simplification_table(
  list(
    LM1.full.0 = LM1.full.0,
    LM1.full.1 = LM1.full.1,
    LM1.full.2 = LM1.full.2,
    LM1.full.3 = LM1.full.3,
    LM1.full.4 = LM1.full.4
  ),
  c(
    "temp_mean_summer:fishing_status",
    "Z_mean:fishing_status",
    "temp_mean_winter:fishing_status",
    "N:fishing_status"
  )
)

## Model 2: phenotypic change
LM2_parameters <- make_parameter_table(LM2.full.5)
LM2_simplification <- make_simplification_table(
  list(
    LM2.full.0 = LM2.full.0,
    LM2.full.1 = LM2.full.1,
    LM2.full.2 = LM2.full.2,
    LM2.full.3 = LM2.full.3,
    LM2.full.4 = LM2.full.4,
    LM2.full.5 = LM2.full.5
  ),
  c(
    "N:fishing_status",
    "temp_mean_summer:fishing_status",
    "temp_mean_winter:fishing_status",
    "N",
    "temp_mean_summer"
  )
)

write.csv(LM1_parameters,
          "table_LM1_population_parameter_estimates.csv",
          row.names = FALSE)
write.csv(LM1_simplification,
          "table_LM1_population_model_simplification.csv",
          row.names = FALSE)
write.csv(LM2_parameters,
          "table_LM2_phenotype_parameter_estimates.csv",
          row.names = FALSE)
write.csv(LM2_simplification,
          "table_LM2_phenotype_model_simplification.csv",
          row.names = FALSE)

#
###

dev.off()
