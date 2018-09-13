
#' Title
#'
#' @param beta
#' @param x.self.design
#' @param z.design
#' @param covar.names
#' @param subtypes.names
#'
#' @return
#' @export
#'
#' @examples
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
