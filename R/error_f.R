

#' Title
#'
#' @param delta_old Input for `delta_old`.
#' @param delta_new Input for `delta_new`.
#'
#' @return A result from this computation.
#' @export
#'
error_f <- function(delta_old,delta_new){
  max(abs(delta_new-delta_old))/(max(abs(delta_old),0.1))
}
