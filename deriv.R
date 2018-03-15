#' coefficient for mixed derivatives
JMcx = c(1,-1,-1,1)
#' coefficient for regular derivatives
JMstep = c(1,-1)

#' compute the mixed derivative of a function
#' 
#' takes two indexes, a function and its argument and compute the mixed derivative 
#' using two 3-points stencils.
#' 
#' @param fn a function that taxes a vector argument and a list of argument
#' @param x the point to compute the mixed derivative
#' @param pos a vector of two indexes for the mixed derivative
#' @param method.args the list of arguments that fn takes
#' @param delta centered difference step value
#' @return mixed derivative of f with x_pos1 and x_pos2
#' @export
JMmixedDeriv = function(fn, x, pos, method.args = list(), delta = 0.01){
  if(1 == length(delta))
    delta = rep(delta, length(x))
  h = 0
  index = 1
  for(p1 in JMstep)
    for(p2 in JMstep){
      xdiff = x
      xdiff[pos[2]] = xdiff[pos[2]] + p2 * delta[pos[2]]
      xdiff[pos[1]] = xdiff[pos[1]] + p1 * delta[pos[1]]
      h = h + JMcx[index] * do.call(fn,c(list(xdiff)
                                         ,method.args))
      index = index + 1
    }
  h = h / (4 * delta[pos[1]] * delta[pos[2]])		
  h
}

#' compute the derivative of a function
#' 
#' takes one index and a function and return df/dx_pos (x) using a 3-point stencil
#' 
#' @param fn a function 
#' @param x the point to estimate the derivative
#' @param pos the index for the derivative
#' @param method.args the arguments of the functions in a list
#' @param delta the centered difference step value
#' @return df/dx_pos (x)
#' @export
JMderiv = function(fn, x, pos, method.args = list(), delta = 0.01){
  if(1 == length(delta))
    delta = rep(delta, length(x))	
  xmin = x
  xmax = x
  xmin[pos] = x[pos] - delta[pos]
  xmax[pos] = x[pos] + delta[pos]
  
  as.numeric((do.call(fn, c(list(xmax),method.args)) 
              - do.call(fn, c(list(xmin),method.args)))/ (2*delta[pos]))
}

#' compute the gradient of a function
#' 
#' compute the gradient of a function with specified centered difference step.
#' 
#' @param fn a function to compute gradient
#' @param x point at which we want the gradient
#' @param method.args the list of arguments the function takes
#' @param delta centered difference step value
#' @return grad f(x)
#' @export
JMgrad = function(fn, x, method.args = list(), delta = 0.01){
  n = length(x)
  g = rep(0,n)
  for(i in 1:n)
    g[i] = JMderiv(fn,x,i,method.args,delta)
  g
}

#' compute the hessian matrix of a function
#' 
#' compute the hessian matrix of a function with specified centered difference step.
#' 
#' @param fn a function to compute hessian matrix
#' @param x point at which we want the hessian matrix
#' @param method.args the list of arguments the function takes
#' @param delta centered difference step value
#' @return H f(x)
#' @export
JMhessian = function(fn, x, method.args = list(), delta = 0.01){
  n = length(x)
  H = matrix(0,n,n)
  if(1 == length(delta))
  	delta = rep(delta, n)
  for(i in 1:n)
    for(j in 1:i)					
      H[i,j] = H[j,i] = JMmixedDeriv(fn,x,c(i,j),method.args,delta)
  H
}
