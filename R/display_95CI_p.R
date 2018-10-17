




#' DisplaySecondStageTestResult
#'
#' @param logodds the log odds ratio vector
#' @param sigma the covariance matrix of the log odds ratio vector
#' @param self.design self design matrix
#' @param places numerical places to keep
#'
#' @keywords internal
#'
DisplaySecondStageTestResult = function(logodds,sigma,self.design=F,
                                        places= 5){
  var.logodds <- diag(sigma)
  logodds.low <- logodds-1.96*sqrt(var.logodds)
  logodds.high <- logodds+1.96*sqrt(var.logodds)
  odds <- exp(logodds)
  odds.low <- exp(logodds.low)
  odds.high <- exp(logodds.high)
  odds <- exp(logodds)
  odds.low <- exp(logodds.low)
  odds.high <- exp(logodds.high)
  places <- 5
  odds <- round(odds,places)
  odds.low <- round(odds.low,places)
  odds.high <- round(odds.high,places)
  p.global.assoc <- GlobalTestForAssoc(logodds,sigma)
  p.global.heter <- GlobalTestForHeter(logodds,sigma,self.design)
  p.individual.heter <- IndividualHeterTest(logodds,sigma)
  result = data.frame(matrix(0,1,2*length(odds)+2))
  for(i in 1:length(logodds)){
    result[1,2*i-1] <- paste0(odds[i],"(",odds.low[i],"-",
                             odds.high[i],")")
    result[1,2*i] <- p.individual.heter[i]
  }
  result[,2*length(odds)+1] <- p.global.assoc
  result[,2*length(odds)+2] <- p.global.heter
  return(result)
}



#' Display the first stage parameters
#'
#' @param logodds the log odds ratio vectors
#' @param sigma the covariance matrix for the log odds ratio vectors
#' @param places the numerical places to keep
#'
#' @keywords internal
#'
DisplayFirstStageTestResult = function(logodds,sigma,
                                       places = 5){
  var.logodds <- diag(sigma)
  logodds.low <- logodds-1.96*sqrt(var.logodds)
  logodds.high <- logodds+1.96*sqrt(var.logodds)
  odds <- exp(logodds)
  odds.low <- exp(logodds.low)
  odds.high <- exp(logodds.high)
  odds <- exp(logodds)
  odds.low <- exp(logodds.low)
  odds.high <- exp(logodds.high)

  odds <- round(odds,places)
  odds.low <- round(odds.low,places)
  odds.high <- round(odds.high,places)
  p.individual.heter <- IndividualHeterTest(logodds,sigma)
  result = data.frame(matrix(0,1,2*length(odds)))
  for(i in 1:length(logodds)){
    result[1,2*i-1] <- paste0(odds[i],"(",odds.low[i],"-",
                              odds.high[i],")")
    result[1,2*i] <- p.individual.heter[i]
  }

  return(result)
}





#' Calcualte the fixed effect score test p value based on score and information matrix
#'
#' @param score the score vector
#' @param infor the information matrix for the score
#' @param places the numerical places to keep
#'

#' @keywords internal
#'

DisplayFixedScoreTestResult <- function(score,infor,
                                        places=5){
  p.value.GTA <- ScoreGlobalTestForAssoc(score,infor)
  p.value.GTA <- round(p.value.GTA*10^(floor(-log10(p.value.GTA))+places))/(10^(floor(-log10(p.value.GTA))+places))
  return(p.value.GTA)
}

#' Calculate the mixed effect model p value based on score and information matrix
#'
#' @param score.fix fixed effect score
#' @param infor.fix fixed effect information matrix
#' @param score.random random effect score
#' @param infor.random random effect information matrix
#'
#' @keywords internal
#'
DisplayMixedScoreTestResult <- function(score.fix,infor.fix,score.random,infor.random){
  p.value.GTH.mixed <- ScoreMixedGlobalTestForHeter(score.random,infor.random)
  p.value.GTA.mixed <- ScoreMixedGlobalTestForAssoc(p.value.GTH.mixed,
                                                    score.fix,
                                                    infor.fix)
  result <- c(p.value.GTA.mixed,p.value.GTH.mixed)
  return(result)
}



#' Calculate the individual test p value
#'
#' @param logodds the log odds ratio vector
#' @param sigma the covariance matrix of the log odds ratio
#' @param places the numerical places to keep
#'
#' @keywords internal
#'
DisplayIndTestResult = function(logodds,sigma,
                                places=5){
  var.logodds <- diag(sigma)
  logodds.low <- logodds-1.96*sqrt(var.logodds)
  logodds.high <- logodds+1.96*sqrt(var.logodds)
  odds <- exp(logodds)
  odds.low <- exp(logodds.low)
  odds.high <- exp(logodds.high)
  places <- 5
  odds <- round(odds,places)
  odds.low <- round(odds.low,places)
  odds.high <- round(odds.high,places)
  p.individual.heter <- IndividualHeterTest(logodds,sigma)
  result = NULL
  for(i in 1:length(logodds)){
    temp <- c(odds,odds.low,odds.high,p.individual.heter[i])
    result= rbind(result,temp)
  }
  return(result)
}
