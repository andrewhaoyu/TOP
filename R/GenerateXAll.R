#' Title
#'
#' @param y
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#'
#' @return
#' @export
#'
#' @examples
GenerateXAll <- function(y,baselineonly,additive,pairwise.interaction,saturated){
  n = nrow(y)
  ###initial x.all to use cbind
  x.all = rep(1,n)
  if(is.null(baselineonly)==0){
    x.all = cbind(x.all,baselineonly)
  }
  if(is.null(additive)==0){
    x.all = cbind(x.all,additive)
  }
  if(is.null(pairwise.interaction)==0){
    x.all = cbind(x.all,pairwise.interaction)
  }
  if(is.null(saturated)==0){
    x.all = cbind(x.all,saturated)
  }
  x.all <- x.all[,-1]
  return(as.matrix(x.all))
}

