
#' Two stage model with self design second stage matrix
#'
#' @param y the phenotype file. The first column is the case control disease status. The other columns are the tumor characteristics status
#' @param x.self.design the covariate you want to use self designed second stage matrix
#' @param z.design self designed second stage matrix
#' @param baselineonly the covariates to be adjusted used baseline effect only model. This assumes the odds ratio of the covariates for all the subtpes to be the same. 
#' @param additive the covariates to be adjusted used the additive two-stage model
#' @param pairwise.interaction the covariates to be adjusted used the pairwise interaction two-stage model
#' @param saturated the covariates to be adjusted used the saturated two-stage model. This model assumes every subtype has their specific odds ratio. It's equivalent to the polytmous model. 
#' @param missingTumorIndicator The indicators to show the tumor characteristics are missing. In the example, we put missing tumor characteristics as 888. Note, for all the controls subjects, they don't have tumor characteristics. So their tumor characteristics are put as NA instead of 888 to differentiate with cases missing tumor characteristics.

#' @param z.all if you want to have differnt self designed second stage matrix for different covariantes, then you can directly construct the second stage matrix for all of the covariates.
#' @param delta0 the starting value for the second stage parameters. By defualt, we will use the empirical distribution of the subtypes.
#'
#'
#' @return the result is a list containing 9 elements. 1. the second stage parameters 2. the covariance matrix for the second stage parameters. 3. the second stage parameters organzied for the self desinged covariate 4. The odds ratio of the self designed subtypes5. Global association test and global heterogeneity test result (Wald test based) 6. The first stage parameter organized for self designed covariates 7. First stage odds test results of all the subtypes. 8. Likelihood 9. AIC
#' @export
#'
#' @examples 
#'#load in the breast cancer example
#' data(data, package="TOP") #load in the breast cancer example
#'#this is a simulated breast cancer example
#'#there are around 5000 breast cancer cases and 5000 controls, i.e. people without disease
#' data[1:5,]
#' 
#'#four different tumor characteristics were included, ER (positive vs negative), PR (positive vs negative), HER2 (positive vs negative), grade (ordinal 1, 2, 3)
#'#the phenotype file
#' y <- data[,1:5]
#'#generate the combinations of all the subtypes
#'#by default, we remove all the subtypes with less than 10 cases
#'z.standard <- GenerateZstandard(y)
#'M <- nrow(z.standard) #M is the total number of first stage subtypes

#'#initial a z.design matrix with M rows, and 5 columns
#'#each row represent a first stage subtype
#'#each column represent an aggregated subtype
#'z.design <- matrix(0,M,5)
#'#define names for the five intrinsic subtypes
#'colnames(z.design) <- c("HR+_HER2-_lowgrade",
#'                        "HR+_HER2+",
#'                        "HR+_HER2-_highgrade",
#'                        "HR-_HER2+", 
#'                        "HR-_HER2-")
#'#To construct a self design second stage matrix,
#'#we need to find the correpsonding first stage subtypes
#'#belonging to specific aggregated subtypes
#'#for first subtype HR+_HER2-_lowgrade
#'idx.1 <- which((z.standard[,1]==1|z.standard[,2]==1)
#'               &z.standard[,3]==0
#'               &(z.standard[,4]==1|z.standard[,4]==2))
#'z.design[idx.1,1] <- 1
#'#for second subtype HR+_HER2+
#'idx.2 <- which((z.standard[,1]==1|z.standard[,2]==1)
#'               &z.standard[,3]==1)
#'z.design[idx.2,2] <- 1
#'#for third subtype HR+_HER2-_highgrade
#'idx.3 <- which((z.standard[,1]==1|z.standard[,2]==1)
#'               &z.standard[,3]==0
#'               &z.standard[,4]==3)
#'z.design[idx.3,3] <- 1
#'#for third subtype HR-_HER2+
#'idx.4 <- which(z.standard[,1]==0&z.standard[,2]==0
#'               &z.standard[,3]==1)
#'z.design[idx.4,4] <- 1
#'#for third subtype HR-_HER2-
#'idx.5 <- which(z.standard[,1]==0&z.standard[,2]==0
#'               &z.standard[,3]==0)
#'z.design[idx.5,5] <- 1
#'#one SNP and one Principal components (PC1) are the covariates
#'SNP <- data[,6,drop=F]
#'PC1 <- data[,7,drop=F]

#'model.3 <- EMmvpolySelfDesign(y,
#'                              x.self.design = SNP,
#'                              z.design = z.design,
#'                              additive=PC1,
#'                              missingTumorIndicator = 888)


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
    y.pheno.complete <- y[-missing.data.vec,drop=F]
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
