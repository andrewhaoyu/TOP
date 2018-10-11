#' Title
#'
#' @param snppro
#' @param n
#'
#' @return
#' @export
#'
#' @examples
convert <- function(snppro,n){
  snpvalue <- rep(0,n)
  temp <- .C("convert",as.integer(n),as.numeric(snppro),as.numeric(snpvalue))

  return(temp[[3]])
}
