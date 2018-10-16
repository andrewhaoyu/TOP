#' Check whether control subtypes tumor characteristics is NA
#'
#' @param y.case.control the case control vector, should be coded as 0, 1. 
#' @param y.tumor the tumor characteristics matrix of all the subjects
#'
#' @return checking result
#' @keywords internal
#'
CheckControlTumor <- function(y.case.control,y.tumor){
  idx <- which(y.case.control==0)
  y.tumor.control = y.tumor[idx,]
  return(any(!is.na(y.tumor.control)))


}
