#' Title
#'
#' @param num
#' @param size
#' @param ind
#'
#' @return
#' @export
#'
#' @examples
startend <- function(num,size,ind){
  split.all <- split(1:num,cut(1:num,size))
  temp <- split.all[[ind]]
  start <- temp[1]
  end <- temp[length(temp)]
  return(c(start,end))
}
