
#' Title
#'
#' @param beta
#' @param covar.names
#' @param subtypes.names
#'
#' @return
#' @export
#'
#' @examples
GenerateFirstStageMat <- function(beta,
                                  covar.names,
                                  subtypes.names
                                  ){

  covar.names <- c("Intercept",covar.names)
  covar.number <- length(covar.names)
  first.mat <- matrix(beta,ncol = covar.number,byrow=T)
  colnames(first.mat) <- covar.names
  rownames(first.mat) <- subtypes.names


  return(first.mat)

}
