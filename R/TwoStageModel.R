

#' Two-stage polytomous logistic regression
#'
#' @param y the phenotype file. The first column is the case control disease status. The other columns are the tumor characteristics status
#' @param baselineonly the covariates to be adjusted used baseline effect only model. This assumes the odds ratio of the covariates for all the subtpes to be the same. 
#' @param additive the covariates to be adjusted used the additive two-stage model
#' @param pairwise.interaction the covariates to be adjusted used the pairwise interaction two-stage model
#' @param saturated the covariates to be adjusted used the saturated two-stage model. This model assumes every subtype has their specific odds ratio. It's equivalent to the polytmous model. 
#' @param missingTumorIndicator The indicators to show the tumor characteristics are missing. In the example, we put missing tumor characteristics as 888. Note, for all the controls subjects, they don't have tumor characteristics. So their tumor characteristics are put as NA instead of 888 to differentiate with cases missing tumor characteristics.
#' @param delta0 the starting value for the second stage parameters. By defualt, we will use the empirical distribution of the subtypes.
#'
#' @return the result is a list containing 9 elements. 1. the second stage parameters 2. the covariance matrix for the second stage parameters. 3. the second stage parameters organzied for each covariate. 4. The case-control odds ratio and case-case odds ratios of tumor characteristics. 5. Global association test and global heterogeneity test result (Wald test based) 6. The first stage parameter organized for each covariate 7. First stage odds ratio of all the subtypes. 8. Likelihood 9. AIC
#' @export
#'
#' @examples
#' 
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
#'#fit the additive two-stage model
#' model.1 <- TwoStageModel(y=y,additive=cbind(SNP,PC1),
#'                          missingTumorIndicator = 888)
#'                          
#' #the second stage odds ratio (95% CI) and p-value. 
#' #the baseline effect is the case-control odds ratio of the reference subtype (ER-,PR-,HER2-,grade0). 
#' #The main effect are the case-case odds ratio of the tumor characteristics.
#' model.1[[4]] 
#' 
#' #the global association test and global heterogeneity test result of the covariates.
#' model.1[[5]] 
#' 
#' model.1[[7]] #the case-control odds ratios for all of the subtypes.
#' 
#'#instead of additive model, you can also try different combinations.
#'# For example, for the PC1, we use the additive model, 
#'#but for SNP, we use the baseline only model.
#' model.2 <- TwoStageModel(y=y,baselineonly = SNP,
#'                       additive=PC1,
#'                       missingTumorIndicator = 888)

TwoStageModel <- function(y,
                          baselineonly=NULL,
                          additive=NULL,
                          pairwise.interaction=NULL,
                          saturated=NULL,
                          missingTumorIndicator = NULL,
                          delta0 = NULL){
  if(is.null(missingTumorIndicator)==1){
    return(Mvpoly(y,
                  baselineonly,
                  additive,
                  pairwise.interaction,
                  saturated,
                  delta0 = delta0))
  }else{

      return(EMmvpoly(y,
                      baselineonly,
                      additive,
                      pairwise.interaction,
                      saturated,
                      missingTumorIndicator,
                      delta0 = delta0))


  }
}
