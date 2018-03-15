#' minimum probability to assign when estimated probability or mass is 0
minProba = .Machine$double.xmin

#' log-likelihood function calculator
#' 
#' calculate the loglikelihood function at point b from a function
#' that gives the probability (or mass) for each observation. corrects for 0 values.
#' 
#' @param b the point we evalute the LL function at
#' @param f the function that calculates the vector of probability (or mass)
#' @param args arguments of the function f in a list
#' @return LL evaluated at b, as specified by function f
#' @export
LLWrapper = function(b, f, args){
  p = f(b,args)
  #print(p)
  p[p < minProba] = minProba
  p[is.nan(p)] = minProba
  LL = sum(log(p))
  #print(LL)
  LL
}

#' calculate the gradient of log-likelihood function
#' 
#' calculate the gradient of LL function at point b from a function
#' that gives the probability (or mass) for each observation. corrects for 0 values.
#' 
#' @param b the point we evalute the gradient at
#' @param f the function that calculates the vector of probability (or mass)
#' @param args arguments of the function f in a list
#' @return grad LL evaluated at b, as specified by function f
#' @export
LLGradWrapper = function(b, f, args){
  if(args[["verbose"]])
    print(LLWrapper(b,f,args))
  delta = args[["delta"]]
  JMgrad(fn = LLWrapper, b, method.args = c(list(f = f), list(args = args)), delta)
}

#' calculate the hessian matrix of log-likelihood function
#' 
#' calculate the gradient of LL function at point b from a function
#' that gives the probability (or mass) for each observation. corrects for 0 values.
#' 
#' @param betaHat the point we evalute the gradient at
#' @param f the function that calculates the vector of probability (or mass)
#' @param args arguments of the function f in a list
#' @return grad LL evaluated at b, as specified by function f
#' @export
LLHessianWrapper = function(betaHat, f, args){
  JMhessian(LLWrapper, betaHat, c(list(f = f), list(args = args))
            , args[["delta"]])
}
