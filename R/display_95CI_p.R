




#' Title
#'
#' @param logodds
#' @param sigma
#' @param self.design
#' @param places
#'
#' @return
#' @export
#'
#' @examples
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



#' Title
#'
#' @param logodds
#' @param sigma
#' @param places
#'
#' @return
#' @export
#'
#' @examples
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





#' Title
#'
#' @param score
#' @param infor
#' @param places
#'
#' @return
#' @export
#'
#' @examples
DisplayFixedScoreTestResult <- function(score,infor,
                                        places=5){
  p.value.GTA <- ScoreGlobalTestForAssoc(score,infor)
  p.value.GTA <- round(p.value.GTA*10^(floor(-log10(p.value.GTA))+places))/(10^(floor(-log10(p.value.GTA))+places))
  return(p.value.GTA)
}

#' Title
#'
#' @param score.baseline
#' @param infor.baseline
#' @param score.casecase
#' @param infor.casecase
#'
#' @return
#' @export
#'
#' @examples
DisplayMixedScoreTestResult <- function(score.baseline,infor.baseline,score.casecase,infor.casecase){
  p.value.GTH.mixed <- ScoreMixedGlobalTestForHeter(score.casecase,infor.casecase)
  p.value.GTA.mixed <- ScoreMixedGlobalTestForAssoc(p.value.GTH.mixed,
                                                    score.baseline,
                                                    infor.baseline)
  result <- c(p.value.GTA.mixed,p.value.GTH.mixed)
  return(result)
}



#' Title
#'
#' @param logodds
#' @param sigma
#' @param places
#'
#' @return
#' @export
#'
#' @examples
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
