

#' Title
#'
#' @param logodds Input for `logodds`.
#' @param sigma Input for `sigma`.
#' @param M Input for `M`.
#' @param second.stage.mat Input for `second.stage.mat`.
#'
#' @return A result from this computation.
#' @export
#'

SecondStageTest <- function(logodds,sigma,M,second.stage.mat){
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

  all.covar.names <- colnames(second.stage.mat)
  all.second.stage.names <- row.names(second.stage.mat)

  covar.names <- NULL
  second.stage.effect.names <- NULL
  for(i in 1:ncol(second.stage.mat)){
    for(j in 1:nrow(second.stage.mat)){
      if(is.na(second.stage.mat[j,i])==FALSE){
        covar.names = c(covar.names,all.covar.names[i])
        second.stage.effect.names =
          c(second.stage.effect.names,all.second.stage.names[j])
      }
    }
  }

  second.n <- length(covar.names)
  odds <- odds[1:second.n]
  odds.low <- odds.low[1:second.n]
  odds.high <- odds.high[1:second.n]
  p.individual.heter <- p.individual.heter[1:second.n]

  result <- data.frame(CovarName = covar.names,
                       SecondStageEffect = second.stage.effect.names,
                       odds,
                       odds.low,
                       odds.high,
                       p.individual.heter)

  colnames(result) <- c("Covariate","SecondStageEffect",
                        "OddsRatio",
                        "OddsRatio(95%CI low)",
                        "OddsRatio(95%CI high)",
                        "Pvalue")
  return(result)

}
