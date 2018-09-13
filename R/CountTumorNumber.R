##count variable number
#' Title
#'
#' @param covar
#'
#' @return
#' @export
#'
#' @examples
CountCovarNumber <- function(covar) {
  if(is.null(covar)){
    return(0)
  }else if(is.vector(covar)){
    return(1)
  }else{
    return(ncol(covar))
  }
}
