
#' Title
#'
#' @param logodds1
#' @param sigma1
#' @param logodds2
#' @param sigma2
#'
#' @return
#' @export
#'
#' @examples
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
#' @param score1
#' @param infor1
#' @param score2
#' @param infor2
#'
#' @return
#' @export
#'
#' @examples
ScoreMetaAnalysis <- function(score1,infor1,score2,infor2){
  infor.meta <- infor1+infor2
  score.meta <- score1+score2

  return(list(score.meta = score.meta,
              infor.meta = infor.meta))

}
