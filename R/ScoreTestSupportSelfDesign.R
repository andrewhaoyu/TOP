###Finish implementing ScoreTestSupport C function into the try5.c code


#' Title
#'
#' @param y
#' @param x.self.design
#' @param z.design
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param missingTumorIndicator
#'
#' @return
#' @export
#'
#' @examples
ScoreTestSupportSelfDesign <- function(y,
                            x.self.design,
                            z.design,
                             baselineonly=NULL,
                             additive=NULL,
                             pairwise.interaction=NULL,
                             saturated=NULL,
                             missingTumorIndicator = 888){
  tumor.number <- ncol(y)-1
  y.case.control <- y[,1]
  y.tumor <- y[,2:(tumor.number+1)]
  y.pheno.complete <- GenerateCompleteYPheno(y,missingTumorIndicator)
  freq.subtypes <- GenerateFreqTable(y.pheno.complete)
  if(CheckControlTumor(y.case.control,y.tumor)==1){
    return(print("ERROR:The tumor characteristics for control subtypes should put as NA"))
  }
  tumor.names <- colnames(y.tumor)
  if(is.null(tumor.names)){
    tumor.names <- paste0(c(1:tumor.number))
  }
  tumor.character.cat = GenerateTumorCharacterCat(y.pheno.complete)
  z.design.baselineonly <- GenerateZDesignBaselineonly(tumor.character.cat,
                                                       tumor.number,
                                                       tumor.names,
                                                       freq.subtypes)
  z.design.additive <- GenerateZDesignAdditive(tumor.character.cat,
                                               tumor.number,
                                               tumor.names,
                                               freq.subtypes)
  z.design.pairwise.interaction <- GenerateZDesignPairwiseInteraction(tumor.character.cat,
                                                                      tumor.number,
                                                                      tumor.names,
                                                                      freq.subtypes)
  z.design.saturated <- GenerateZDesignSaturated(tumor.character.cat,
                                                 tumor.number,
                                                 tumor.names,
                                                 freq.subtypes)
  full.second.stage.names <- colnames(z.design)
  covar.names <- GenerateSelfCovarName(x.self.design,
                                       baselineonly,
                                       additive,
                                       pairwise.interaction,
                                       saturated)
  z.all <- ZSelfDesigntoZall(x.self.design,
                             baselineonly,
                             additive,
                             pairwise.interaction,
                             saturated,
                             z.design,
                             z.design.baselineonly,
                             z.design.additive,
                             z.design.pairwise.interaction,
                             z.design.saturated)
  delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
  #x.all has no intercept yet
  #we will add the intercept in C code
  x.all <- GenerateSelfXAll(y,x.self.design,baselineonly,additive,pairwise.interaction,saturated)
  ###z standard matrix means the additive model z design matrix without baseline effect
  ###z standard matrix is used to match the missing tumor characteristics to the complete subtypes

  y <- as.matrix(y)
  x.all <- as.matrix(x.all)
  z.standard <- z.design.additive[,-1]
  M <- as.integer(nrow(z.standard))
  p.main <- ncol(z.standard)+1

  Score.Support = EMStepScoreTestSupport(delta0,y,x.all,z.standard,z.all,missingTumorIndicator)


  # score_support_result <- score_support(pxx,x.all,baselineonly,z.all,z.standard,y_em)
  #score_test_mis <- score_test_mis(y_em,baselineonly,score_support_result)
  #return(list(score_c=score_test_mis$score_c,infor_c = score_test_mis$infor_c))
  return(Score.Support)

}
