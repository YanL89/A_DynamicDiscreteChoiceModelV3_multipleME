#' translate a dynamic specification into a logit specification with data
#' 
#' @param dynSpec the dynamic specification
#' @export
#' @return the spec and the data
dyn2logitV2 = function(dynSpec){
  generic = dynSpec$generic
  specific = dynSpec$specific
  
  #prepend a "zero" because the "keep" option is a full fledged alternative now
  if(0 != length(generic))
		for(i in 1:length(generic))
			generic[[i]] = c("zero",generic[[i]])
  
  keepvars = c(dynSpec$payoff_alt, dynSpec$payoff_time, dynSpec$payoff_global)
  specific = c(list(keepvars), specific)
  
  dynSpec = checkFillDynSpec(dynSpec)
  D = dyn2logitV2data(dynSpec) 
  
  specLogit = list(
    generic = generic,
    specific = specific,
    Y = "logitY",
    ASC = dynSpec$ASC)
  
  list(specLogit = specLogit, D = D)
}

#' compute the data for the dynamic to logit conversion
#' 
#' @param spec the dynamic spec
#' @export
#' @return a unified dataset
dyn2logitV2data = function(spec){
  # first time period
  Z = spec$First
  Dt = spec$modifyD(spec$D, spec$Time, spec$Global, Z, 1)
  Z = updateZ(Z,spec$D[[1]],spec$C[,1],spec$varNames)
  D = Dt
  for(t in 2:spec$nTime){
    Dt = spec$modifyD(spec$D, spec$Time, spec$Global, Z, t)
    Z = updateZ(Z,spec$D[[t]],spec$C[,t],spec$varNames)
    D = rbind(D,Dt)
  }
  D$zero = 0
  D$one = 1
  D$logitY = dyn2logitV2choices(spec$C)
  D
}

#' adapts the choice matrix for a logit
#' 
#' @param C the matrix of choices
#' @return a vector of logit choices
dyn2logitV2choices = function(C){
  C2 = C[,1]
  for(i in 2:ncol(C))
    C2 = c(C2, C[,i])
  C2
}
