#' dyn model functions with AR(1) simulation
dynVAR1 = list(
  LLVec = function(b, args){
    Xsim = args$Xsim
    nparH = ncol(args$M[[1]][[1]])
    
    alpha = b[1:nparH]
    beta = b[(nparH+1):length(b)]
    
    LL = 0
    B = args$spec$B
    
    H = list()
    U = list()
    V = list()
    PAllObs = rep(1,args$spec$nObs)

    for(t in 1:args$spec$nTime){
    
    Wtally = rep(0,args$spec$nObs)
    Rtally = rep(0,args$spec$nObs)
    U1tally = matrix(0,args$spec$nObs,args$spec$nAlt)
    
    # simulation of dynamic attributes here
		row = seq(from = 1, to = args$spec$nObs*args$spec$nAlt, by = args$spec$nAlt)
		nSimVar = length(args$colIndexX)   #number of simulated variable
		
    for(bsim in 1:B){
			
      for(l in 1:(args$spec$nLook+1)){
        H[[l]] = as.matrix(args$M[[t]][[l]]) %*% alpha
        if(0 == ncol(H[[l]]))
        	H[[l]] = rep(0, args$spec$nObs)
        
        # fill Xl with the l-th replication of EACH dynamic variable	
        Xl = args$X[[t]][[l]]
        for(v in 1:length(args$colIndexX))
					Xl[args$rowIndexX[v] - 1 + row, args$colIndexX[v]] = Xsim[[t]][[bsim]][[v]][,l]
					
        U[[l]] = matrix(Xl %*% beta, args$spec$nObs, args$spec$nAlt, byrow = T)
        V[[l]] = UtoV(U[[l]])
      }
      # tally W, U1 and R
      Wtally = Wtally + getW(H,V)
      Rtally = Rtally + getR(U[[1]])
      U1tally = U1tally + U[[1]]
      
      } # END OF SIMULATION
      
      W = Wtally / B
      R = Rtally / B
      U1 = U1tally / B
    
      P = getDynProbas(W,R,U1,args$C[,t])
      P = maxVec(P,args$stop[,t])
      PAllObs = PAllObs * P
    }
    #for-loop and print(LL)
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
    
    Zlist = list()
    Zlist[[1]] = Z
    for(t in 1:spec$nTime){
      M[[t]] = list()
      X[[t]] = list()
      for(l in 1:(spec$nLook+1)){
        M[[t]][[l]] = getM(Z, t+l-1, spec)
        Dt = spec$modifyD(spec$D, spec$Time, spec$Global, Z, t+l-1)
        X[[t]][[l]] = create_X(spec$generic, spec$specific, Dt)
      }
      Z = updateZ(Z,spec$D[[t]],spec$C[,t],spec$varNames)
      Zlist[[t+1]] = Z
    }
    
    oneX = X[[1]][[1]]
    # dyn variable index
    colIndexX = c()
    rowIndexX = c()
    for(i in spec$dynvar){
      colIndexX = c(colIndexX, which(colnames(oneX) == i))
      #Yan add 
      for (j in 1:spec$nAlt){
        if(i %in% spec$specific[[j]])
          rowIndexX = c(rowIndexX,j)   
      }
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
    
    list(spec = spec, M = M, X = X, C = spec$C, stop = stopMat, Zlist = Zlist, colIndexX = colIndexX, rowIndexX = rowIndexX)
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
    set.seed(1234) 
    nTime = args$spec$nTime
    B = args$spec$B
    nLook = args$spec$nLook
    nSimVar = length(args$colIndexX)
    nObs = args$spec$nObs
    ARpars = args$spec$dynpars[[1]]
    
    row = seq(from = 1, to = args$spec$nObs*args$spec$nAlt, by = args$spec$nAlt)
    
    Xsim = list()
    Xv = list()
    for(t in 1:nTime){
      Xsim[[t]] = list()
      for(bsim in 1:B){
        Xsim[[t]][[bsim]] = list()
        for(v in 1:nSimVar){
          Xv[[v]] = matrix(0,nObs,nLook+1)
          Xv[[v]][,1] = args$X[[t]][[1]][row + args$rowIndexX[v] - 1,args$colIndexX[v]]
        }
        
        for(i in 1:nObs){
          for(j in 2:(nLook+1)){
            x = Xv[[1]][i,j-1]
            for(v in 2:nSimVar){
              x = cbind(x, Xv[[v]][i,j-1]) 
            }
            y = dynSimVAR(x,ARpars[1:nSimVar],ARpars[(nSimVar+1):(nSimVar*(1+nSimVar))],ARpars[(nSimVar*(1+nSimVar)+1): (nSimVar*(1+2*nSimVar))])
            for(v in 1:nSimVar){
            Xv[[v]][i,j] = y[v] 
            }
          }
        }
        
        for(v in 1:nSimVar){
          Xsim[[t]][[bsim]][[v]] = Xv[[v]]
        }
    }
    }
    list(Xsim = Xsim)
}
) # dynAR1 function list
