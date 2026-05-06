###############
## FUNCTIONS ##
###############

multigrep = function(x, patterns)
{
  s = NULL
  for (pattern in patterns)
  {
    s_ = grep(x = x, pattern = pattern)  
    s = c(s, s_)
  }
  return(s)
}

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