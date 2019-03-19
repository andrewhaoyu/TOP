
#' Generate the support list for score test
#' To construct a score test, we need to first fit the model under the null hypothesis. This function fitts the model under the null hypothesis, then saved all the necessary results into a list to pass to the score test function.
#' @param y the phenotype file. The first column is the case control disease status. The other columns are the tumor characteristics status
#' @param baselineonly the covariates to be adjusted used baseline effect only model. This assumes the odds ratio of the covariates for all the subtpes to be the same. 
#' @param additive the covariates to be adjusted used the additive two-stage model
#' @param pairwise.interaction the covariates to be adjusted used the pairwise interaction two-stage model
#' @param saturated the covariates to be adjusted used the saturated two-stage model. This model assumes every subtype has their specific odds ratio. It's equivalent to the polytmous model. 
#' @param missingTumorIndicator The indicators to show the tumor characteristics are missing. In the example, we put missing tumor characteristics as 888. Note, for all the controls subjects, they don't have tumor characteristics. So their tumor characteristics are put as NA instead of 888 to differentiate with cases missing tumor characteristics.
#'
#' @return return a list for score test function
#' @export
#'
#' @examples
#' data(data, package="TOP") #load in the breast cancer example
#'#this is a simulated breast cancer example
#'#there are around 5000 breast cancer cases and 5000 controls, i.e. people without disease
#' data[1:5,]
#' 
#'#four different tumor characteristics were included, ER (positive vs negative), PR (positive vs negative), HER2 (positive vs negative), grade (ordinal 1, 2, 3)
#'#the phenotype file
#' y <- data[,1:5]
#' 
#'#one SNP and one Principal components (PC1) are the covariates
#' SNP <- data[,6,drop=F]
#' PC1 <- data[,7,drop=F]
#' 
#'#fit the additive two-stage model under the null hypothesis that the second stage parameters of SNP is 0
#'score.support <- ScoreTestSupportMixedModel(y=y,
#'                additive=PC1,
#'                missingTumorIndicato#'r=888)

ScoreTestSupportMixedModel <- function(y,
                             baselineonly=NULL,
                             additive=NULL,
                             pairwise.interaction=NULL,
                             saturated=NULL,
                             missingTumorIndicator = 888,
                             delta0 = NULL){

  y <- as.matrix(y)
  tumor.number <- ncol(y)-1
  y.case.control <- y[,1,drop=F]
  y.tumor <- y[,2:(tumor.number+1),drop=F]
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
  
  if(tumor.number>=2){
    z.design.pairwise.interaction <- GenerateZDesignPairwiseInteraction(tumor.character.cat,
                                                                        tumor.number,
                                                                        tumor.names,
                                                                        freq.subtypes)
    z.design.saturated <- GenerateZDesignSaturated(tumor.character.cat,
                                                   tumor.number,
                                                   tumor.names,
                                                   freq.subtypes)
    
  }else{
    z.design.pairwise.interaction <- z.design.additive
    z.design.saturated <- z.design.additive
    
  }
  z.all <- ZDesigntoZall(baselineonly,
                         additive,
                         pairwise.interaction,
                         saturated,
                         z.design.baselineonly,
                         z.design.additive,
                         z.design.pairwise.interaction,
                         z.design.saturated)

  if(is.null(delta0)==T){
    delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
  }else{
    delta0 =delta0
  }

  #x.all has no intercept yet
  #we will add the intercept in C code
  x.all <- GenerateXAll(y,baselineonly,additive,pairwise.interaction,saturated)
  ###z standard matrix means the additive model z design matrix without baseline effect
  ###z standard matrix is used to match the missing tumor characteristics to the complete subtypes

  z.standard <- z.design.additive[,-1,drop=F]

  Score.Support = EMStepScoreTestSupportMixedModel(delta0,y,x.all,z.standard,z.all,missingTumorIndicator)


  # score_support_result <- score_support(pxx,x.all,baselineonly,z.all,z.standard,y_em)
  #score_test_mis <- score_test_mis(y_em,baselineonly,score_support_result)
  #return(list(score_c=score_test_mis$score_c,infor_c = score_test_mis$infor_c))
  result <- Score.Support
  result[[6]] <- z.design.baselineonly
  result[[7]] <- z.design.additive
  result[[8]] <- z.design.pairwise.interaction
  result[[9]] <- z.design.saturated
  result[[10]] <- z.standard
  return(result)

}
















