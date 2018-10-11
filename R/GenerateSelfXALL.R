
#' Title
#'
#' @param y
#' @param x.self.design
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#'
#' @return
#' @export
#'
#' @examples
GenerateSelfXAll <- function(y,x.self.design,baselineonly,additive,pairwise.interaction,saturated){
  n = nrow(y)
  ###initial x.all to use cbind
  x.all = rep(1,n)
  if(is.null(x.self.design)==0){
    x.all <- cbind(x.all,x.self.design)
  }
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
  x.all <- x.all[,-1,drop=F]
  return(as.matrix(x.all))
}

