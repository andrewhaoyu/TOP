

#' Title
#'
#' @param delta_old
#' @param delta_new
#'
#' @return
#' @export
#'
#' @examples
error_f <- function(delta_old,delta_new){
  max(abs(delta_new-delta_old))/(max(abs(delta_old),0.1))
}
