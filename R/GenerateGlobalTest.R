#' Title
#'
#' @param delta
#' @param sigma
#' @param M
#' @param second.stage.mat
#'
#' @return
#' @export
#'
#' @examples
GenerateGlobalTest <- function(delta,sigma,M,second.stage.mat){
  ind.start=1
  ind.end = 0

  p <- ncol(second.stage.mat)
  covar.names <- colnames(second.stage.mat)
  global.test.for.assoc <- rep(0,p)
  global.test.for.heter <- rep(0,p)
  for(i in 1:p){
    for(j in 1:nrow(second.stage.mat)){
      if(is.na(second.stage.mat[j,i])==F){
        ind.end <- ind.end+1
      }
    }
      delta.covar.i <- delta[ind.start:ind.end]
      sigma.covar.i <- sigma[ind.start:ind.end,
                             ind.start:ind.end]
      global.test.for.assoc[i] <- GlobalTestForAssoc(delta.covar.i,sigma.covar.i)
      global.test.for.heter[i] <- GlobalTestForHeter(delta.covar.i,sigma.covar.i)
      GlobalTestForHeter(delta.covar.i,sigma.covar.i)
      ind.start = ind.end+1

  }

    result <- data.frame(covar.names,
                         global.test.for.assoc,
                         global.test.for.heter
    )
    colnames(result) <- c("Covariate",
                          "Global test for association",
                          "Global test for heterogeneity")
    return(result)

  }
