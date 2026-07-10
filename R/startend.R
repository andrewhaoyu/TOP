#' Title
#'
#' @param num Input for `num`.
#' @param size Input for `size`.
#' @param ind Input for `ind`.
#'
#' @return A result from this computation.
#' @export
#'
startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}
