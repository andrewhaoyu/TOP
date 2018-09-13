##check whether control subtypes tumor characteristics is NA
#' Title
#'
#' @param y.case.control
#' @param y.tumor
#'
#' @return
#' @export
#'
#' @examples
CheckControlTumor <- function(y.case.control,y.tumor){
  idx <- which(y.case.control==0)
  y.tumor.control = y.tumor[idx,]
  return(any(!is.na(y.tumor.control)))


}
