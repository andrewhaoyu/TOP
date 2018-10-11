
#' Title
#'
#' @param y_em
#' @param p
#' @param nparm
#'
#' @return
#' @export
#'
#' @examples
LogLikelihoodwithAIC <- function(y_em,p,nparm){
  n <- nrow(y_em)
  M <- ncol(y_em)



  loglikelihood <- ComputeLogLikelihood(y_em,p)
  AIC <- 2*nparm-2*loglikelihood
  return(list(loglikelihood=loglikelihood,


              AIC = AIC
            ))
}

#' Title
#'
#' @param y_em
#' @param p
#'
#' @return
#' @export
#'
#' @examples
ComputeLogLikelihood <- function(y_em,p){

  n <- nrow(y_em)
  M <- ncol(y_em)
  p.o <- rep(0,n)
  y_em_ind <- t(apply(y_em,1, function(x){
    x[x>1e-09] <- 1
    return(x)
  }))

  loglikelihood <- 0
  for(i in 1:n){
    p.temp <- p[i+n*(0:(M-1))]
    y.temp <- y_em_ind[i,]
    if(any(y.temp==1)==T){
      p.o[i] <- crossprod(y.temp,(p.temp))
    }else{
      p.o[i] <- (1-sum(p.temp))
    }
    loglikelihood <- sum(log(p.o))

    #loglikelihood <- loglikelihood+
    #crossprod(y.temp,log(p.temp))+(1-sum(y.temp))*log(1-sum(p.temp))
  }
  return(loglikelihood)
}


