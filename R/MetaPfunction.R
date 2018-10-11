#' Title
#'
#' @param icog_onco_score_infor_one
#' @param second.num
#'
#' @return
#' @export
#'
#' @examples
MetaPfunction <- function(icog_onco_score_infor_one,second.num){
  score.icog <- as.numeric(icog_onco_score_infor_one[1:(second.num)])
  infor.icog <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1):
                                                              (second.num+second.num^2) ]),
                       ncol = second.num)
  start <- second.num+second.num^2
  score.onco <- as.numeric(icog_onco_score_infor_one[(1+start):
                                                       (second.num+start)])
  infor.onco <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1+start):
                                                              (second.num+second.num^2+start) ]),ncol=second.num)


  meta.result <- ScoreMetaAnalysis(score.icog,infor.icog,
                                   score.onco,infor.onco)
  score.meta <- t(meta.result[[1]])
  infor.meta <- meta.result[[2]]
  DisplayFixedScoreTestResult(score.meta,infor.meta)
}



#' Title
#'
#' @param icog_onco_score_infor_one
#' @param icog_onco_score_infor_casecase_one
#' @param fixed.second.num
#' @param random.second.num
#'
#' @return
#' @export
#'
#' @examples
MetaMixedPfunction <- function(icog_onco_score_infor_one,icog_onco_score_infor_casecase_one,fixed.second.num,random.second.num){
  score.icog <- as.numeric(icog_onco_score_infor_one[1:(fixed.second.num)])
  infor.icog <- matrix(as.numeric(icog_onco_score_infor_one[(fixed.second.num+1):
                                                              (fixed.second.num+fixed.second.num^2) ]),
                       ncol = fixed.second.num)
  start <- fixed.second.num+fixed.second.num^2
  score.onco <- as.numeric(icog_onco_score_infor_one[(1+start):
                                                       (fixed.second.num+start)])
  infor.onco <- matrix(as.numeric(icog_onco_score_infor_one[(fixed.second.num+1+start):
                                                              (fixed.second.num+fixed.second.num^2+start) ]),ncol=fixed.second.num)

  meta.result.fixed <- ScoreMetaAnalysis(score.icog,infor.icog,
                                         score.onco,infor.onco)


  score.meta.fixed <- t(meta.result.fixed[[1]])
  infor.meta.fixed <- meta.result.fixed[[2]]







  score.icog.casecase <- as.numeric(icog_onco_score_infor_casecase_one[1:(random.second.num)])
  infor.icog <- matrix(as.numeric(icog_onco_score_infor_casecase_one[(random.second.num+1):
                                                              (random.second.num+random.second.num^2) ]),
                       ncol = random.second.num)
  start <- random.second.num+random.second.num^2
  score.onco <- as.numeric(icog_onco_score_infor_casecase_one[(1+start):
                                                                (random.second.num+start)])
  infor.onco <- matrix(as.numeric(icog_onco_score_infor_casecase_one[(random.second.num+1+start):
                                                                       (random.second.num+random.second.num^2+start) ]),ncol=random.second.num)

  meta.result.random <- ScoreMetaAnalysis(score.icog,infor.icog,
                                          score.onco,infor.onco)
  score.meta.random <- t(meta.result.random[[1]])
  infor.meta.random <- meta.result.random[[2]]


  result <-   DisplayMixedScoreTestResult(score.meta.fixed,
                                          infor.meta.fixed,
                                          score.meta.random,
                                          infor.meta.random)
  return(result[1])
}











#' Title
#'
#' @param icog_onco_score_infor_one
#' @param icog_onco_score_infor_casecase_one
#' @param fixed.second.num
#' @param random.second.num
#'
#' @return
#' @export
#'
#' @examples
MetaMixedPfunction_temp <- function(icog_onco_score_infor_one,icog_onco_score_infor_casecase_one,fixed.second.num,random.second.num){
  score.icog <- as.numeric(icog_onco_score_infor_one[1:(fixed.second.num)])
  infor.icog <- matrix(as.numeric(icog_onco_score_infor_one[(fixed.second.num+1):
                                                              (fixed.second.num+fixed.second.num^2) ]),
                       ncol = fixed.second.num)
  start <- fixed.second.num+fixed.second.num^2
  score.onco <- as.numeric(icog_onco_score_infor_one[(1+start):
                                                       (fixed.second.num+start)])
  infor.onco <- matrix(as.numeric(icog_onco_score_infor_one[(fixed.second.num+1+start):
                                                              (fixed.second.num+fixed.second.num^2+start) ]),ncol=fixed.second.num)

  meta.result.fixed <- ScoreMetaAnalysis(score.icog,infor.icog,
                                         score.onco,infor.onco)


  score.meta.fixed <- t(meta.result.fixed[[1]])
  infor.meta.fixed <- meta.result.fixed[[2]]




  score.icog <- rep(0,random.second.num)
  #score.icog <- rep(0,temp.n)

  infor.icog <- matrix(0,nrow= random.second.num,
                       ncol = random.second.num)
  start <- random.second.num+random.second.num^2
  score.onco <- as.numeric(icog_onco_score_infor_casecase_one[(1+start):
                                                                (random.second.num+start)])
  infor.onco <- matrix(as.numeric(icog_onco_score_infor_casecase_one[(random.second.num+1+start):
                                                                       (random.second.num+random.second.num^2+start) ]),ncol=random.second.num)

  meta.result.random <- ScoreMetaAnalysis(score.icog,infor.icog,
                                          score.onco,infor.onco)
  score.meta.random <- t(meta.result.random[[1]])
  infor.meta.random <- meta.result.random[[2]]


result <-   DisplayMixedScoreTestResult(score.meta.fixed,
                              infor.meta.fixed,
                              score.meta.random,
                              infor.meta.random)
return(result[1])
}




#' Title
#'
#' @param icog_onco_score_infor_one
#' @param second.num
#'
#' @return
#' @export
#'
#' @examples
MetaFixedPfunction_temp <- function(icog_onco_score_infor_one,second.num){
  log.icog <- as.numeric(icog_onco_score_infor_one[1:(second.num)])
  sigma.icog <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1):
                                                              (second.num+second.num^2) ]),
                       ncol = second.num)
  start <- second.num+second.num^2
  log.onco <- as.numeric(icog_onco_score_infor_one[(1+start):
                                                       (second.num+start)])
  sigma.onco <- matrix(as.numeric(icog_onco_score_infor_one[(second.num+1+start):
                                                              (second.num+second.num^2+start) ]),ncol=second.num)

  meta.result.fixed <- LogoddsMetaAnalysis(log.icog,sigma.icog,
                                         log.onco,sigma.onco)


  logodds.meta.fixed <- meta.result.fixed[[1]]
  sigma.meta.fixed <- meta.result.fixed[[2]]




  p.value <- DisplaySecondStageTestResult(logodds.meta.fixed,sigma.meta.fixed,self.design=T,
                                          places= 5)

  return(list(logodds.meta.fixed,sigma.meta.fixed,
              p.value))
}
