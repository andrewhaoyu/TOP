#' Title
#'
#' @param delta Input for `delta`.
#' @param covariance.delta Input for `covariance.delta`.
#' @param M Input for `M`.
#' @param tumor.names Input for `tumor.names`.
#' @param z.all Input for `z.all`.
#' @param covar.names Input for `covar.names`.
#'
#' @return A result from this computation.
#' @export
#'
TakeoutIntercept <- function(delta,covariance.delta,
                             M,
                             tumor.names,
                             z.all,covar.names){
  beta <- z.all%*%delta
  covariance.beta <-
    z.all%*%covariance.delta%*%t(z.all)
  delta.no.inter <- delta[(M+1):length(delta)]
  covariance.delta.no.inter <-
    covariance.delta[(M+1):length(delta),
                     (M+1):length(delta)]
  p.no.inter <- length(covar.names)
  p <- p.no.inter+1
  beta.no.inter <- beta[-(1+p*(0:(M-1)))]
  covariance.beta.no.inter <- covariance.beta[-(1+p*(0:(M-1))),
                                              -(1+p*(0:(M-1)))]
  result <- list(beta = beta,
                 covariance.beta = covariance.beta,
                 delta.no.inter = delta.no.inter,
                 covariance.delta.no.inter = covariance.delta.no.inter,
                 beta.no.inter = beta.no.inter,
                 covariance.beta.no.inter = covariance.beta.no.inter
                 )
  return(result)

}
