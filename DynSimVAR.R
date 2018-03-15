library("mvtnorm")
#' dynamic variable simulator
#' 
#' use stochastic diffusion process such as AR(1) to simulate dynamic varibles
#' 
#' @param x value of dynamic variable at time t-1
#' @param a slope
#' @param e constant
#' @param s standard deviation of white noise
#' @return value of dynamic variable at time t

dynSimAR = function(x,a,e,s){
  n = length(x)
  a*x + e + rnorm(n,0,s)
}

dynSimVAR = function(x,alpha,betain,s){
  beta= matrix(betain, nrow=length(alpha))
  beta= t(beta)
  sigma = matrix(s, nrow=length(alpha))
  n = nrow(x)
  simX = matrix(0, nrow=n, ncol=length(alpha))
  mean = rep(0, length(alpha))
  for(i in 1:n)
    simX[i,] = alpha + beta%*%x[i,]+t(rmvnorm(n=1, mean=mean, sigma=sigma))
  simX
}