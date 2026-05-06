######
##  ##
######

## Goal:

## Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

###################
## FUNCTIONS SLP ##
###################

## goal: functions to define a single layer perceptron

## activation functions ##
lin       = function(x) x
ddx.lin   = function(x) 1
pol2      = function(x) x^2
ddx.pol2  = function(x) 2*x
sigmo     = function(x) 1/(1+exp(-x)) 
ddx.sigmo = function(x) sigmo(x) * (1 - sigmo(x))
expo      = function(x) exp(x)
ddx.expo  = function(x) exp(x)
relu      = function(x) (x>0)*x
ddx.relu  = function(x) (x>0)
tanh      = function(x) (exp(2*x)-1)/(exp(2*x)+1) 
ddx.tanh  = function(x) (2*exp(2*x))/(exp(2*x)+1) + (exp(2*x)-1)*(-1)*(2*exp(2*x))/((exp(2*x)+1)^2)

## SLP ##
## goal: single layer perceptron returning single output with input x and parameter vector Omega 
# x       - vector - input variables
# Omega   - vector - parameters
# f_sigma - func - activation function
SLP = function(x,Omega,f_sigma) 
{    
    Omega = matrix(Omega,ncol=2 + length(x))
    return(t(Omega[,1])%*%f_sigma(Omega[,2] + Omega[,-c(1:2)]%*%t(t(x))))
}

## ddOmega.SLP ##
## goal: compute the derivative of single layer perceptron wtr to each parameter 
# x           - vector - input variables
# Omega       - vector - parameters
# f_sigma     - func - activation function
# ddu.f_sigma - func - derivative of activation function
ddOmega.SLP = function(x,Omega,f_sigma,ddu.f_sigma)
{
    Omega = matrix(Omega,ncol=2 + length(x))
    x = t(x)
    Omega_1 = Omega[,1]
    Omega_2 = Omega[,2]
    Omega_3 = Omega[,-c(1:2)]
    ddOmega_1 = f_sigma(Omega_2 + Omega_3%*%t(x))
    ddOmega_2 = Omega_1 * ddu.f_sigma(Omega_2 + Omega_3%*%t(x))
    ddOmega_3 = Omega_1%*%x * as.vector(ddu.f_sigma(Omega_2 + Omega_3%*%t(x)))
    return(c(ddOmega_1,ddOmega_2,ddOmega_3))
}

## ddx.SLP ##
## goal: compute the derivative of single layer perceptron wtr to each input variable
# x           - vector - input variables
# Omega       - vector - parameters
# f_sigma     - func - activation function
# ddu.f_sigma - func - derivative of activation function
ddx.SLP = function(x,Omega,f_sigma,ddu.f_sigma)
{
    Omega = matrix(Omega,ncol=2 + length(x))
    x = t(x)
    Omega_1 = Omega[,1]
    Omega_2 = Omega[,2]
    Omega_3 = Omega[,-c(1:2)]
    ddx = Omega_1%*%(Omega_3*as.vector(ddu.f_sigma(Omega_2 + Omega_3%*%t(x))))
    return(ddx)
}

#
###