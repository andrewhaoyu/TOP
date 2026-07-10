logit_inver <- function(x){
  return(exp(x)/(1+exp(x)))
}

.top_message <- function(...) {
  if (isTRUE(getOption("TOP.verbose", FALSE))) {
    message(...)
  }
}

load_so = function(so_file) {
  system.file(so_file, package = "TOP")
}


#' @useDynLib TOP, .registration = TRUE, .fixes = "C_"
#' @importFrom stats pchisq pnorm
#' @importFrom utils capture.output combn object.size
NULL
