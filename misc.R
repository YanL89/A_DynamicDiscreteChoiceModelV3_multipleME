#' some big number for bounds
bigNum = 999999999

#' print the member of a list to the screen
#'
#' @param l a list
#' @param name the name of an element
#' @export
printIfElem = function(l,name){
  if(name %in% names(l))
    cat(name, " : ", l[[name]], "\n")
}

#' fetches an element in a list if it exists
#' 
#' returns an element from a list or a default value if the element does not exist
#' 
#' @param l the list
#' @param name the name of the element we want
#' @param default the default value if the element does not exist
#' @return the element or default
#' @export
getElem = function(l, name, default){
  if(name %in% names(l))
    return(l[[name]])
  return(default)
}

#' vector wise logical and
#' 
#' returns true if all elements of a vector are 1
#' 
#' @param b the vector of booleans
#' @return true if all elements of b are true, false otherwise
#' @export
logicalAnd = function(b){
  length(b) == sum(b)
}

#' get starting values for model
#'
#' @param madelFns madel function specifications
#' @param spec model specification
#' @param D the data
#' @return the starting value for the solver
#' @export
getStart = function(modelFns, spec, D){
  start = NULL
  if("start" %in% names(spec)){
    start = scan(spec[["start"]])
  } else {
    start = modelFns$computeStart(spec, D)
  }
  start  
} 

#' get lower bound for solver
#' 
#' @param misc a list that may contain LB
#' @param start the starting value
#' @return LB = 0 or LB specified in misc
#' @export
getLb = function(misc, start){
  lb = rep(-bigNum, length(start))
  if("lb" %in% names(misc))
    lb = misc[["lb"]]
  lb
}

#' get upper bound for solver
#' 
#' @param misc a list that may contain UB
#' @param start the starting value
#' @return UB = 0 or UB specified in misc
#' @export
getUb = function(misc, start){
  ub = rep(-bigNum, length(start))
  if("ub" %in% names(misc))
    ub = misc[["ub"]]
  ub
}

#' compute standard deviations for model parameters
#' 
#' @param spec model specifications
#' @param betaHat the estimated parameters
#' @param modelFns list of model functions
#' @param args model arguments for LL and stuff
#' @return SD of betaHat
#' @export
getSD = function(spec, betaHat, modelFns, args){
  SD = 0
  if(spec[["SD"]] == "none")
    SD = NoneSD(betaHat)
  if(spec[["SD"]] == "hessian")
    SD = HessianSD(betaHat, modelFns, args)
  if(spec[["SD"]] == "bootstrap")
    SD = BootstrapSD(betaHat, args, modelFns, spec)#, start, lb, ub)
  SD
}
