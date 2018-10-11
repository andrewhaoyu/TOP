
#' Title
#'
#' @param y.pheno.complete
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param x.self.design
#' @param z.design
#'
#' @return
#' @export
#'
#' @examples
#'
#'
InitialSetup <- function(y.pheno.complete,
                         baselineonly,
                         additive,
                         pairwise.interaction,
                         saturated,
                         x.self.design= NULL,
                         z.design = NULL
                         ){
  if(is.null(x.self.design)==T){
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
    z.design.pairwise.interaction <- GenerateZDesignPairwiseInteraction(tumor.character.cat,
                                                                        tumor.number,
                                                                        tumor.names,
                                                                        freq.subtypes)
    z.design.saturated <- GenerateZDesignSaturated(tumor.character.cat,
                                                   tumor.number,
                                                   tumor.names,
                                                   freq.subtypes)
    full.second.stage.names <- colnames(z.design.saturated)
    covar.names <- GenerateCovarName(baselineonly,
                                     additive,
                                     pairwise.interaction,
                                     saturated,
                                     x.self.design = x.self.design)

    z.all <- ZDesigntoZall(baselineonly,
                           additive,
                           pairwise.interaction,
                           saturated,
                           z.design.baselineonly,
                           z.design.additive,
                           z.design.pairwise.interaction,
                           z.design.saturated)
    z.standard <- z.design.additive[,-1]
    delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
    return(list(delta0 = delta0,z.all=z.all,z.standard= z.standard,z.deisign.baselineonly = z.design.baselineonly,z.design.additive=z.design.additive,z.design.pairwise.interaction=z.design.pairwise.interaction,z.design.saturated=z.design.saturated,covar.names = covar.names,tumor.names=tumor.names
    ))
  }else{
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
    z.design.pairwise.interaction <- GenerateZDesignPairwiseInteraction(tumor.character.cat,
                                                                        tumor.number,
                                                                        tumor.names,
                                                                        freq.subtypes)
    z.design.saturated <- GenerateZDesignSaturated(tumor.character.cat,
                                                   tumor.number,
                                                   tumor.names,
                                                   freq.subtypes)
    full.second.stage.names <- colnames(z.design.saturated)
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
    z.standard <- z.design.additive[,-1]
    delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
    return(list(delta0 = delta0,z.all=z.all,z.standard= z.standard,z.deisign.baselineonly = z.design.baselineonly,z.design.additive=z.design.additive,z.design.pairwise.interaction=z.design.pairwise.interaction,z.design.saturated=z.design.saturated,covar.names = covar.names,tumor.names=tumor.names
    ))
  }

}







#' Title
#'
#' @param model.result
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param z.standard
#' @param covar.names
#' @param delta
#' @param z.design.additive
#' @param z.design.pairwise.interaction
#' @param z.design.saturated
#' @param tumor.names
#' @param z.all
#' @param x.self.design
#' @param z.design
#'
#' @return
#' @export
#'
#' @examples
SummaryResult <- function(model.result,
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
                          z.all,
                          x.self.design = NULL,
                          z.design = NULL
                          ){

if(is.null(x.self.design)){
  full.second.stage.names <- colnames(z.design.saturated)
  M <- as.integer(nrow(z.standard))
  ###delta represent second stage parameters
  delta <- model.result$delta
  covariance.delta <- solve(model.result$infor_obs)
  loglikelihood <- model.result$loglikelihood
  AIC <- model.result$AIC
  second.stage.mat <-
    GenerateSecondStageMat(baselineonly,
                           additive,
                           pairwise.interaction,
                           saturated,
                           M,
                           full.second.stage.names,
                           covar.names,
                           delta,
                           z.design.additive,
                           z.design.pairwise.interaction,
                           z.design.saturated)
  ##take out the intercept from second stage parameters

  takeout.intercept.result <- TakeoutIntercept(delta,covariance.delta,
                                               M,
                                               tumor.names,
                                               z.all,covar.names)
  beta <- takeout.intercept.result$beta
  covariance.beta <- takeout.intercept.result$covariance.beta
  delta.no.inter <- takeout.intercept.result$delta.no.inter
  covariance.delta.no.inter <-
    takeout.intercept.result$covariance.delta.no.inter
  beta.no.inter <- takeout.intercept.result$beta.no.inter
  covariance.beta.no.inter <- takeout.intercept.result$covariance.beta.no.inter



  second.stage.test <- SecondStageTest(delta.no.inter,covariance.delta.no.inter,M,second.stage.mat)
  global.test <- GenerateGlobalTest(delta.no.inter,
                                    covariance.delta.no.inter,
                                    M,
                                    second.stage.mat)
  ##beta represent first stage parameters

  subtypes.names <- GenerateSubtypesName(z.design.additive,M,
                                         tumor.names)
  first.stage.mat <- GenerateFirstStageMat(beta,
                                           covar.names,
                                           subtypes.names)

  first.stage.test <- FirstStageTest(beta.no.inter,
                                     covariance.beta.no.inter,
                                     M,
                                     first.stage.mat)
  return(list(delta=delta,covariance.delta=covariance.delta,second.stage.mat = second.stage.mat,second.stage.test,global.test,first.stage.mat,first.stage.test,loglikelihood = loglikelihood,
              AIC = AIC,beta=beta,covariance.beta=covariance.beta,
              z.standard=z.standard))
}else{
  if(is.null(colnames(z.design))){
    full.second.stage.names <- paste0("self_design_group",c(1:ncol(z.design)))
  }else{
    full.second.stage.names <- colnames(z.design)
  }

  M <- as.integer(nrow(z.standard))
  ###delta represent second stage parameters
  delta <- model.result$delta
  covariance.delta <- solve(model.result$infor_obs)
  loglikelihood <- model.result$loglikelihood
  AIC <- model.result$AIC
  second.stage.mat <-
    GenerateSelfSecondStageMat(x.self.design,
                               z.design,
                               M,
                               full.second.stage.names,
                               delta)
  ##take out the intercept from second stage parameters

  takeout.intercept.result <- TakeoutIntercept(delta,covariance.delta,
                                               M,
                                               tumor.names,
                                               z.all,
                                               covar.names)
  beta <- takeout.intercept.result$beta
  covariance.beta <- takeout.intercept.result$covariance.beta
  delta.no.inter <- takeout.intercept.result$delta.no.inter
  covariance.delta.no.inter <-
    takeout.intercept.result$covariance.delta.no.inter
  beta.no.inter <- takeout.intercept.result$beta.no.inter
  covariance.beta.no.inter <- takeout.intercept.result$covariance.beta.no.inter



  second.stage.test <- SecondStageTest(delta.no.inter,covariance.delta.no.inter,M,second.stage.mat)
  global.test <- GenerateGlobalTest(delta.no.inter,
                                    covariance.delta.no.inter,
                                    M,
                                    second.stage.mat)
  ##beta represent first stage parameters

  subtypes.names <- GenerateSubtypesName(z.design.additive,M,
                                         tumor.names)
  first.stage.mat <- GenerateSelfFirstStageMat(beta,
                                               x.self.design,
                                               z.design,
                                               covar.names,
                                               subtypes.names
  )
  first.stage.test <- SelfFirstStageTest(beta.no.inter,
                                         covariance.beta.no.inter,
                                         M,
                                         first.stage.mat,
                                         covar.names)

  return(list(delta=delta,covariance.delta=covariance.delta,second.stage.mat = second.stage.mat,second.stage.test,global.test,first.stage.mat,first.stage.test,loglikelihood = loglikelihood,
              AIC = AIC,beta=beta,covariance.beta=covariance.beta,
              z.standard=z.standard))
}

}



