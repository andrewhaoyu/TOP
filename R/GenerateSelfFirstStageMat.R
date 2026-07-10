
#' Title
#'
#' @param beta Input for `beta`.
#' @param x.self.design Input for `x.self.design`.
#' @param z.design Input for `z.design`.
#' @param covar.names Input for `covar.names`.
#' @param subtypes.names Input for `subtypes.names`.
#'
#' @return A result from this computation.
#' @export
#'
GenerateSelfFirstStageMat <- function(beta,
                                  x.self.design,
                                  z.design,
                                  covar.names,
                                  subtypes.names
){
  self.design.number <- CountCovarNumber(x.self.design)
  covar.names <- c("Intercept",covar.names)
  covar.number <- length(covar.names)
  first.mat <- matrix(beta,ncol = covar.number,byrow=T)
  colnames(first.mat) <- covar.names
  rownames(first.mat) <- subtypes.names
  first.mat <- first.mat[,c(2:(1+self.design.number)),drop=F]

  return(first.mat)

}
