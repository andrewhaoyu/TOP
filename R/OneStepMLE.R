
#' Title
#'
#' @param y
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
OneStepMLE <- function(y,
                     baselineonly,
                     additive,
                     pairwise.interaction,
                     saturated,
                     missingTumorIndicator){
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

  delta0 = initial.set$delta0
  z.all = initial.set$z.all
  z.standard = initial.set$z.standard
  z.deisign.baselineonly = initial.set$z.design.baseline.only
  z.design.additive = initial.set$z.design.additive
  z.design.pairwise.interaction = initial.set$z.design.pairwise.interaction
  z.design.saturated = initial.set$z.design.saturated
  x.all <- as.matrix(GenerateXAll(y,baselineonly,additive,pairwise.interaction,saturated))
  covar.names <- initial.set$covar.names


  prob.fit.result <- ProbFitting(delta0,y.pheno.complete,x.all.complete,z.standard,z.all)
  y.fit <- prob.fit.result[[1]]
  M <- as.integer(nrow(z.standard))
  p.main <- ncol(z.standard)+1

  tol <- as.numeric(1e-04)


  delta_old <- delta0

  N <- as.integer(nrow(x.all.complete))

  M <- as.integer(nrow(z.standard))

  NCOV   <- as.integer(ncol(x.all.complete))
  NM     <- N*M
  nparm  <- as.integer(length(delta0))
  deltai <- as.numeric(delta0)

  NITER  <- as.integer(500)
  Y <- as.numeric(as.vector(y.fit.complete))
  X <- as.numeric(as.vector(x.all.complete))
  ZallVec = as.numeric(as.vector(z.all))
  Znr = as.integer(nrow(z.all))
  Znc = as.integer(ncol(z.all))
  debug     <- as.integer(1)
  ret_rc    <- as.integer(1)
  ret_delta <- as.numeric(rep(-9999, nparm))
  ret_info <- as.numeric(rep(-9999,nparm^2))
  ret_p <- as.numeric(rep(0,NM))
  ret_lxx <- as.numeric(rep(0,NM))
  loglikelihood <- as.numeric(-1);



  temp <- .C("Mvpoly_complete",deltai, nparm, Y=Y, X, ZallVec,Znr,Znc, N, M, NCOV, NITER, tol,
             debug, ret_rc=ret_rc, ret_delta=ret_delta,ret_info=ret_info,ret_p=ret_p,loglikelihood = loglikelihood)


  delta0 <- temp$ret_delta


  y.pheno.misonly <-   y[missing.data.vec,]
  x.all.misonly <- x.all[missing.data.vec,]
  prob.fit.result <- ProbFitting(delta0,y,x.all,z.standard,z.all)
  y.fit <- prob.fit.result[[1]]

  missing.vec <- as.numeric(as.vector(prob.fit.result[[2]]))
  missing.mat <- prob.fit.result[[3]]
  missing.mat.vec <- as.numeric(as.vector(missing.mat))
  missing.number <- as.integer(length(missing.vec))
  complete.vec <- prob.fit.result[[4]]
  y.fit.complete <- y.fit[complete.vec,]





  N <- as.integer(nrow(x.all))
  NM     <- N*M

  deltai <- as.numeric(delta0)



  Y <- as.numeric(as.vector(y.fit))
  X <- as.numeric(as.vector(x.all))
  ret_p <- as.numeric(rep(0,NM))
  ret_lxx <- as.numeric(rep(0,NM))



  temp <- .C("OneStepMLE",deltai, nparm, Y=Y, X, ZallVec,Znr,Znc, N, M, NCOV, NITER, tol,
             debug, ret_rc=ret_rc, ret_delta=ret_delta,ret_info=ret_info,ret_p=ret_p,missing.vec,
             missing.mat.vec,missing.number,loglikelihood = loglikelihood)


  info <- matrix(unlist(temp$ret_info),nparm,nparm)
  result <- list(temp$ret_delta,info,
                 temp$ret_p)
  y_em <- matrix(unlist(temp$Y),N,M)

  # infor_mis_c <- infor_mis(y_em,x.all,z.all)
  #infor_obs <- result[[2]]-infor_mis_c
  delta=result[[1]]
  infor_obs=result[[2]]
  p=result[[3]]
  loglikelihood = temp$loglikelihood
  AIC = 2*nparm - 2*loglikelihood

  OneStep.result <- list(delta=delta,
                         infor_obs=infor_obs,
                         p=p,y_em=y_em,
                         M=M,
                         NumberofTumor=ncol(z.standard),
                         loglikelihood = loglikelihood,
                         AIC = AIC
  )

  covariance.delta <- solve(OneStep.result$infor_obs)
  loglikelihood <- OneStep.result$loglikelihood
  AIC <- OneStep.result$AIC
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

}
