#' Complete tho argument lists if some elements are not there
#' 
#' @param args the list of arguments
#' @param newParam the new parameter to add if not present
#' @param defVal the default value
#' @export
addVal = function(args, newParam, defVal){
	if( ! newParam %in% names(args))
		args[[newParam]] = defVal
  args
}
