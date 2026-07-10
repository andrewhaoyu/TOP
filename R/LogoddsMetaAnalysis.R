
#' Title
#'
#' @param logodds1 Input for `logodds1`.
#' @param sigma1 Input for `sigma1`.
#' @param logodds2 Input for `logodds2`.
#' @param sigma2 Input for `sigma2`.
#'
#' @return A result from this computation.
#' @export
#'
LogoddsMetaAnalysis <- function(logodds1,sigma1,logodds2,sigma2){
  sigma1.inv <- solve(sigma1)
  sigma2.inv <- solve(sigma2)
  sigma.meta <- solve(sigma1.inv+sigma2.inv)
  logodds.meta <- sigma.meta%*%(sigma1.inv%*%logodds1+sigma2.inv%*%logodds2)

  return(list(logodds.meta = logodds.meta,
         sigma.meta = sigma.meta))

}


#' Title
#'
#' @param score1 Input for `score1`.
#' @param infor1 Input for `infor1`.
#' @param score2 Input for `score2`.
#' @param infor2 Input for `infor2`.
#'
#' @return A result from this computation.
#' @export
#'
ScoreMetaAnalysis <- function(score1,infor1,score2,infor2){
  meta <- top_meta(list(
    list(score = score1, information = infor1),
    list(score = score2, information = infor2)
  ))

  return(list(score.meta = meta$score,
              infor.meta = meta$information))

}
