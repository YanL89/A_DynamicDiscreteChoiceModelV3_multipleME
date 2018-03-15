#' return unambiguous values for no standard deviations
#' 
#' return 9999999 values when we do not want to compute standard deviations
#' 
#' @param beta_hat the final estimate
#' @return a vector of 9999999 of the same size
#' @export
NoneSD = function(betaHat){
  rep(9999999, length(betaHat))
}

#' hessian matrix based standard deviations
#' 
#' computes the standard deviations of the model based on the hessian
#' matrix evaluated at the maximum of the LL function.
#' 
#' @param betaHat maximum found
#' @param modelFns list of model functions
#' @param args arguments of the LL functions
#' @return estimated SD of beta hat
#' @export
HessianSD = function(betaHat, modelFns, args){
  H = LLHessianWrapper(betaHat, modelFns$LLVec, args)
  sqrt(diag(solve(-H)))	
}	

#' compute boobstrat LL function with weights
#' 
#' computes the bootstrap LL function with weights corresponding to how many
#' times each specific observation has been chosen in the bootstrap sample
#' 
#' @param b parameter value
#' @param f function that computes the probability vector (components of LL)
#' @param args arguments of the LL function f
#' @param w vector of weights for each observation (integer values)
#' @return bootstrap LL function evaluated at b
#' @export
LLWrapperBoot = function(b, f, args, w){
  p = f(b,args)
  p[p < minProba] = minProba
  # w tell how many times we should count each obs
  sum(w * log(p))
}

#' compute gradient of boobstrat LL function with weights
#' 
#' computes the gradient of the bootstrap LL function with weights corresponding 
#' to how many times each specific observation has been chosen in the bootstrap sample
#' 
#' @param b parameter value
#' @param f function that computes the probability vector (components of LL)
#' @param args arguments of the LL function f
#' @param w vector of weights for each observation (integer values)
#' @return grad (bootstrap LL function) evaluated at b
#' @export
LLGradWrapperBoot = function(b, f, args, w){
  # print(LLWrapper(b,f,args))
  delta = args[["delta"]]
  JMgrad(f = LLWrapperBoot, b, c(list(f = f), list(args = args)
                                 , list(w = w)), delta)
}

#' compute non parametric bootstrap variance estimation of a model
#' 
#' compute bootstrap variance estimation by non-parametric resampling of the sample
#' 
#' @param b the estimated maximum
#' @param args arguments of model functions
#' @param modelFns functions that specify the model, in a list
#' @param spec specification of the model (the thing the user supplies)
#' @return estimated standard deviations
#' @export
BootstrapSD = function(b, args, modelFns, spec){
  n = length(modelFns$LLVec(b,args))
  nboot = spec$nboot
  all_estimates = matrix(0,length(b), nboot)
  for(i in 1:nboot){
    #if(0 == i %% 10)
    cat("bootstrap iteration ",i,"\n")
    bootsamp = sample(n,n,replace=TRUE)
    w = rep(0,n)
    for(j in bootsamp)
      w[j] = w[j] + 1
    start_boot = rep(0, length(b))
    O = optim(fn = LLWrapperBoot, gr = LLGradWrapperBoot
              , method = "BFGS", hessian =F, par = start_boot
              , control = list(fnscale = -1, reltol = spec[["reltol"]])
              , args = args, w = w, f = modelFns$LLVec)
    all_estimates[,i] = O$par
    
  }
  sqrt(diag(var(t(all_estimates))))
}
