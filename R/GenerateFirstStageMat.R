
#' Title
#'
#' @param beta Input for `beta`.
#' @param covar.names Input for `covar.names`.
#' @param subtypes.names Input for `subtypes.names`.
#'
#' @return A result from this computation.
#' @export
#'
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
