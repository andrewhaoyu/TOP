
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
#' @param z.all
#' @param delta0
#'
#'
#' @return
#' @export
#'
#' @examples

EMmvpolySelfDesign <- function(y,
                               x.self.design,
                               z.design,
                               baselineonly=NULL,
                               additive=NULL,
                               pairwise.interaction=NULL,
                               saturated=NULL,
                               missingTumorIndicator = 888,
                               z.all=NULL,
                               delta0 = NULL){
  if(is.null(z.all)){
    missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator)
    y.pheno.complete <- y[-missing.data.vec,]
    initial.set <- InitialSetup(y.pheno.complete,
                                baselineonly,
                                additive,
                                pairwise.interaction,
                                saturated,
                                x.self.design,
                                z.design
    )
    ###z standard matrix means the additive model z design matrix without baseline effect
    ###z standard matrix is used to match the missing tumor characteristics to the complete subtypes
    if(is.null(delta0)){
      delta0 = initial.set$delta0
    }
    z.all = initial.set$z.all
    z.standard = initial.set$z.standard
    z.deisign.baselineonly = initial.set$z.design.baseline.only
    z.design.additive = initial.set$z.design.additive
    z.design.pairwise.interaction = initial.set$z.design.pairwise.interaction
    z.design.saturated = initial.set$z.design.saturated
    x.all <- GenerateSelfXAll(y,x.self.design,baselineonly,additive,pairwise.interaction,saturated)
    covar.names <- initial.set$covar.names
    tumor.names <- initial.set$tumor.names

    ###z standard matrix means the additive model z design matrix without baseline effect
    ###z standard matrix is used to match the missing tumor characteristics to the complete subtypes

    y <- as.matrix(y)
    x.all <- as.matrix(x.all)

    M <- as.integer(nrow(z.standard))
    p.main <- ncol(z.standard)+1

    model.result = EMStep(delta0,as.matrix(y),x.all,z.standard,z.all,missingTumorIndicator)
    ###delta represent second stage parameters
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
                                    z.all,
                                    x.self.design,
                                    z.design
    )

    #   pxx = EM.result[[3]]
    #   y_em = EM.result[[4]]
    #  score_support_result <- score_support(pxx,x.all,baselineonly,z.all,z.standard,y_em)
    #  #return(score_support_result)
    # score_test_mis_result <- score_test_mis(y_em,baselineonly,score_support_result)

    return(summary.result)
  }else{

  }
}
