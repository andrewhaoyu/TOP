logit_inver <- function(x){
  return(exp(x)/(1+exp(x)))
}

load_so = function(so_file) {
  system.file(so_file, package = "TOP")
}


#' @useDynLib TOP
#' @importFrom Rcpp sourceCpp
NULL
