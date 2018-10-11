#' Title
#'
#' @param y
#' @param missingTumorIndicator
#'
#' @return
#' @export
#'
#' @examples
GenerateZstandard <- function(y,
                              missingTumorIndicator = 888){
  if(is.null(missingTumorIndicator)){
  y.pheno.complete <- y
  tumor.number <- ncol(y)-1
  y.case.control <- y[,1]
  y.tumor <- y[,2:(tumor.number+1)]
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
    y.case.control <- y[,1]
    y.tumor <- y[,2:(tumor.number+1)]
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
#' Title
#'
#' @param y
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param missingTumorIndicator
#' @param delta0
#'
#' @return
#' @export
#'
#' @examples
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
