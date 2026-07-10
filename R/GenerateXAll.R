#' Title
#'
#' @param y Input for `y`.
#' @param baselineonly Input for `baselineonly`.
#' @param additive Input for `additive`.
#' @param pairwise.interaction Input for `pairwise.interaction`.
#' @param saturated Input for `saturated`.
#'
#' @return A result from this computation.
#' @export
#'
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

