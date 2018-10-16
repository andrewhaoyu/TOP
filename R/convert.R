#' Convert the impute 2 data into the 0, 1, 2 coding
#'
#' @param snppro impute 2 data with 3 probabilities for every SNP (aa, Aa, AA)
#' @param n the number of subjects
#'

#' @keywords internal
#'
convert <- function(snppro,n){
  snpvalue <- rep(0,n)
  temp <- .C("convert",as.integer(n),as.numeric(snppro),as.numeric(snpvalue))

  return(temp[[3]])
}
