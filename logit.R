#' logit model functions
logit = list(
  LLVec = function(b, args){
    n = length(args$Y)
    nalt = max(args$Y) - min(args$Y) + 1	
    
    expU = exp(args$X %*% b)
    
    num = expU[seq(from = 1, length = n, by = nalt) + args$Y]
    denum = rowSums(matrix(expU, ncol = nalt, byrow = T))
    
    num / denum
  },
  
  computeArgs = function(spec, D){
    X = create_X(spec$generic, spec$specific, D, spec$ASC)
    Y = D[[spec$Y]]
    Y = Y - min(Y)
    delta = getElem(spec, name = "delta", default = 0.001) 
    list(X = X, Y = Y, delta = delta)
  },
  
  computeStart = function(spec, D){
  	numVar = getNumVar(spec$generic, spec$specific, D) + 
  			length(formatASC(spec$ASC, length(spec$specific)))
    rep(0, numVar)
  },
  
  computeMisc = function(spec, D){
    list(names = getNames(spec$generic, spec$specific, D, spec$ASC))
  },
  
  apply = function(spec, Val, DVal){
   X = create_X(spec$common, spec$specific, DVal)
   U = matrix(X %*% Val$results$beta_hat, nrow = nrow(DVal), byrow = T)
   U = exp(U)
   P = U / rowSums(U)
   
   list(expected = colSums(P), observed = table(DVal[,spec$Y]))
  }
) # logit function list

#' generate a logit choice
genLogit = function(generic, specific, D, b){
  nObs = nrow(D)
  nAlts = length(specific)
  X = create_X(generic, specific, D)
  U = matrix(X %*% b, nObs, nAlts, byrow = T)
  U = U + matrix(rgumbel(nObs * nAlts), nObs, nAlts)
  apply(U, 1, which.max)
}

logitApplyP = function(M,spec,D){
  nAlt = length(spec$specific)
  X = create_X(spec$common, spec$specific, D)
  U = matrix( X %*% M$results$beta_hat, ncol = nAlt, byrow = T)
  U = exp(U)
  denom = rowSums(U)
  for(i in 1:ncol(U))
    U[,i] = U[,i] / denom
  U
}

#' computes the logit probas from Utilities
#' 
#' @param U utility matrix
#' @export
#' @return matrix of choice probabilities
getPAlts = function(U){
  U = exp(U)
  U / rowSums(U)
}
