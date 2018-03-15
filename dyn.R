#' dyn model functions
dyn = list(
  LLVec = function(b, args){
    nparH = ncol(args$M[[1]][[1]])
    
    alpha = b[1:nparH]
    beta = b[(nparH+1):length(b)]
    
    LL = 0
    
    H = list()
    U = list()
    V = list()
    PAllObs = rep(1,args$spec$nObs)

    for(t in 1:args$spec$nTime){
      for(l in 1:(args$spec$nLook+1)){
        H[[l]] = as.matrix(args$M[[t]][[l]]) %*% alpha
        if(0 == ncol(H[[l]]))
        	H[[l]] = rep(0, args$spec$nObs)
        U[[l]] = matrix(args$X[[t]][[l]] %*% beta, args$spec$nObs, args$spec$nAlt, byrow = T)
        V[[l]] = UtoV(U[[l]])    #max one from U
      }
      W = getW(H,V)  #obtain expected payoff at current time t
      R = getR(U[[1]])   #sum of log(utilities) of all alternatives at time t, mode of the utility
      P = getDynProbas(W,R,U[[1]],args$C[,t])   #get probability of keeping or choosing certain alternative
      P = maxVec(P,args$stop[,t])  #if out_of_market (alter=1), let the probability = 1
      
      PAllObs = PAllObs * P
    }
   
    for(i in 1:args$spec$nObs){
    LL = LL + log(PAllObs[i])
    }
    #print(LL)
    PAllObs
  },
  
  computeArgs = function(spec, D){
    Z = spec$First
    M = list()
    X = list()
    
    for(t in 1:spec$nTime){
      M[[t]] = list()
      X[[t]] = list()
      for(l in 1:(spec$nLook+1)){
        M[[t]][[l]] = getM(Z, t+l-1, spec)  #get all attributes for pay-off at t
        Dt = spec$modifyD(spec$D, spec$Time, spec$Global, Z, t+l-1)   #get all attributes at t
        X[[t]][[l]] = create_X(spec$generic, spec$specific, Dt)
      }
      Z = updateZ(Z,spec$D[[t]],spec$C[,t],spec$varNames)
    }
    
    # check if we reached stop point already
    stopMat = matrix(0, spec$nObs, spec$nTime)  
    
    if(spec$outTime > 0){
      for(i in 1:spec$nObs){
        t = 2
       while(t <= spec$nTime){ 
          if(spec$C[i,t-1] %in% spec$stopAlt){
            stopMat[i,t:min(t+spec$outTime-1, spec$nTime)] = (spec$C[i,t-1] %in% spec$stopAlt)   
            t = min(t+spec$outTime-1, spec$nTime) + 2
          } else {t = t+1}
         }
        }
    }
    
    list(spec = spec, M = M, X = X, C = spec$C, stop = stopMat)   
  },

  
  computeStart = function(spec, D){
    npar = length(spec$payoff_alt) + length(spec$payoff_time) + length(spec$payoff_global)
    npar = npar + length(spec$generic)
    for(s in spec$specific)
      npar = npar + length(s)
    rep(0, npar)
  },
  
  computeMisc = function(spec, D){
    n = c(spec$payoff_alt, spec$payoff_time, spec$payoff_global)
    comIndex = 1
    for(com in spec$generic){
      n = c(n,paste("common",comIndex,sep="_"))
      comIndex = comIndex + 1
    }
    for(s in spec$specific)
      n = c(n,s)
    list(names = n)
  },
  
  computeSimData = function(args, D){
    list()
  }
) # dyn function list


dynApply = function(M,spec){
  args = dyn$computeArgs(spec,NULL)
  b = M$results$beta_hat
  predCount = matrix(0, spec$nTime, spec$nAlt+1)
  
  nparH = ncol(args$M[[1]][[1]])
  alpha = b[1:nparH]
  beta = b[(nparH+1):length(b)]
  
  H = list()
  U = list()
  V = list()
  for(t in 1:args$spec$nTime){
    for(l in 1:(args$spec$nLook+1)){
      H[[l]] = as.matrix(args$M[[t]][[l]]) %*% alpha
      if(0 == ncol(H[[l]]))
        H[[l]] = rep(0, args$spec$nObs)
      U[[l]] = matrix(args$X[[t]][[l]] %*% beta, args$spec$nObs, args$spec$nAlt, byrow = T)
      V[[l]] = UtoV(U[[l]])
    }
    W = getW(H,V)
    R = getR(U[[1]])
    for(choice in 0:spec$nAlt){
      P = getDynProbas(W,R,U[[1]],rep(choice,spec$nObs))
      P = minVec(P,1 - args$stop[,t])
      predCount[t,choice+1] = sum(P)    
    }
  }
  predCount
}


dynCount = function(spec){
  count = matrix(0, spec$nTime, spec$nAlt + 1)
  for(t in 1:spec$nTime)
    for(a in 0:spec$nAlt)
      count[t,a+1] = sum(spec$C[,t] == a)
  count
}
