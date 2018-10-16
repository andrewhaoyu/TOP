#' Count variable number
#'
#' @param covar the covariate matrix
#'
#' @keywords internal
#'
CountCovarNumber <- function(covar) {
  if(is.null(covar)){
    return(0)
  }else if(is.vector(covar)){
    return(1)
  }else{
    return(ncol(covar))
  }
}
