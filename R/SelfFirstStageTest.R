
#' Title
#'
#' @param beta
#' @param sigma
#' @param M
#' @param first.stage.mat
#'
#' @return
#' @export
#'
#' @examples
SelfFirstStageTest <- function(beta,sigma,M,first.stage.mat,covar.names){
  M <- nrow(first.stage.mat)
  self.design.number <- ncol(first.stage.mat)
  self.covar.names <- colnames(first.stage.mat)

  logodds <- beta
  ###take out intercept
  all.covar.names <- colnames(first.stage.mat)
  all.subtypes.names <- row.names(first.stage.mat)
  first.stage.mat <- first.stage.mat
  if(is.vector(first.stage.mat)){
    first.stage.mat <- matrix(first.stage.mat,ncol=1)
  }
  ind.delta <- 0
  ind.covar <- 0
  var.logodds <- diag(sigma)
  logodds.low <- logodds-1.96*sqrt(var.logodds)
  logodds.high <- logodds+1.96*sqrt(var.logodds)
  odds <- exp(logodds)
  odds.low <- exp(logodds.low)
  odds.high <- exp(logodds.high)
  places <- 2
  odds <- round(odds,places)
  odds.low <- round(odds.low,places)
  odds.high <- round(odds.high,places)
  p.individual.heter <- IndividualHeterTest(logodds,sigma)



  covar.names <- NULL
  subtypes.names <- NULL
  for(i in 1:ncol(first.stage.mat)){
    for(j in 1:nrow(first.stage.mat)){
      if(is.na(first.stage.mat[j,i])==F){
        covar.names = c(covar.names,all.covar.names[i])
        subtypes.names =
          c(subtypes.names,all.subtypes.names[j])
      }
    }
  }


  odds <- odds[1:(M*self.design.number)]
  odds.low <- odds.low[1:(M*self.design.number)]
  odds.high <- odds.high[1:(M*self.design.number)]
  p.individual.heter <- p.individual.heter[1:(M*self.design.number)]


  result <- data.frame(CovarName = self.covar.names,
                       Subtypes = subtypes.names,
                       odds,
                       odds.low,
                       odds.high,
                       p.individual.heter)

  colnames(result) <- c("Covariate","Subtypes",
                        "OddsRatio",
                        "OddsRatio(95%CI low)",
                        "OddsRatio(95%CI high)",
                        "Pvalue")
  return(result)

}

