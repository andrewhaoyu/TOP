#' Generate the combinations of all the tumor characteristics.
#'
#' @param y the phenotype file. The first column is the case control disease status. The other columns are the tumor characteristics status
#' @param missingTumorIndicator The indicators to show the tumor characteristics are missing. In the example, we put missing tumor characteristics as 888. Note, for all the controls subjects, they don't have tumor characteristics. So their tumor characteristics are put as NA instead of 888 to differentiate with cases missing tumor characteristics.
#'
#' @return a matrix with all the combinations of tumor characteristics. All the subtypes with less than 10 are removed by default.
#' @export
#'
#' @examples
#' data(data, package="TOP") #load in the breast cancer example
#'#this is a simulated breast cancer example
#'#there are around 5000 breast cancer cases and 5000 controls, i.e. people without disease
#' data[1:5,]

#'#four different tumor characteristics were included, ER (positive vs negative), PR (positive vs negative), HER2 (positive vs negative), grade (ordinal 1, 2, 3)
#'#the phenotype file
#'y <- data[,1:5]
#'#generate the combinations of all the subtypes
#'#by default, we remove all the subtypes with less than 10 cases
#'z.standard <- GenerateZstandard(y)

GenerateZstandard <- function(y,
                              missingTumorIndicator = 888){
  if(is.null(missingTumorIndicator)){
  y.pheno.complete <- y
  tumor.number <- ncol(y)-1
  y.case.control <- y[,1,drop=F]
  y.tumor <- y[,2:(tumor.number+1),drop=F]
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
  z.standard <- z.design.additive[,-1,drop=F]
  return(z.standard)

  }else{
    missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator)
    y.pheno.complete <- y[-missing.data.vec,]
    y <- y.pheno.complete
    tumor.number <- ncol(y)-1
    y.case.control <- y[,1,drop=F]
    y.tumor <- y[,2:(tumor.number+1),drop=F]
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
    z.standard <- z.design.additive[,-1,drop=F]
    return(z.standard)

  }

}
#' Two-stage model MLE estimation with EM algorithm
#'
#' @param y the phenotype file. The first column is the case control disease status. The other columns are the tumor characteristics status
#' @param baselineonly the covariates to be adjusted used baseline effect only model. This assumes the odds ratio of the covariates for all the subtpes to be the same. 
#' @param additive the covariates to be adjusted used the additive two-stage model
#' @param pairwise.interaction the covariates to be adjusted used the pairwise interaction two-stage model
#' @param saturated the covariates to be adjusted used the saturated two-stage model. This model assumes every subtype has their specific odds ratio. It's equivalent to the polytmous model. 
#' @param missingTumorIndicator The indicators to show the tumor characteristics are missing. In the example, we put missing tumor characteristics as 888. Note, for all the controls subjects, they don't have tumor characteristics. So their tumor characteristics are put as NA instead of 888 to differentiate with cases missing tumor characteristics.
#' @param delta0 the starting value for the second stage parameters. By defualt, we will use the empirical distribution of the subtypes.
#'
#' @keywords internal
#'
EMmvpoly <- function(y,
                          baselineonly=NULL,
                          additive=NULL,
                          pairwise.interaction=NULL,
                          saturated=NULL,
                          missingTumorIndicator = 888,
                     delta0= NULL){

  missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator)
  y.pheno.complete <- y[-missing.data.vec,]
  initial.set <- InitialSetup(y.pheno.complete,
                           baselineonly,
                           additive,
                           pairwise.interaction,
                           saturated
  )
  ###z standard matrix means the additive model z design matrix without baseline effect
  ###z standard matrix is used to match the missing tumor characteristics to the complete subtypes
  if(is.null(delta0)==T){
    delta0 = initial.set$delta0
  }else{
    delta0 = delta0
  }

  z.all = initial.set$z.all
  z.standard = initial.set$z.standard
  z.deisign.baselineonly = initial.set$z.design.baseline.only
  z.design.additive = initial.set$z.design.additive
  z.design.pairwise.interaction = initial.set$z.design.pairwise.interaction
  z.design.saturated = initial.set$z.design.saturated
  x.all <- as.matrix(GenerateXAll(y,baselineonly,additive,pairwise.interaction,saturated))
  covar.names <- initial.set$covar.names
  tumor.names <- initial.set$tumor.names


    model.result = EMStep(delta0,as.matrix(y),x.all,z.standard,z.all,missingTumorIndicator)

    summary.result <- SummaryResult(model.result,
                                    baselineonly,
                                    additive,
                                    pairwise.interaction,
                                    saturated,
                                    z.standard,
                                    covar.names,
                                    delta,
                                    z.design.additive,
                                    z.design.pairwise.interaction,
                                    z.design.saturated,
                                    tumor.names,
                                    z.all
    )


    return(summary.result)

}
