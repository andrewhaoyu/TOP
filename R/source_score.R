fixedEffect <- function(data, cc.var, subtype.vars, test.vars, covars.obj=NULL, 
                 options=NULL) {

  options      <- check_options(options, fixed=TRUE)
  subtype.vars <- unique(subtype.vars)
  test.vars    <- unique(test.vars)
  checkDataAndVars(data, cc.var, subtype.vars, test.vars, covars.obj)

  # Get all objects for no missing values in test.vars
  obj0   <- getAllObjects(data, cc.var, subtype.vars, test.vars, covars.obj, options)
  subset <- obj0$subset

  ret <- list()
  for (i in 1:length(test.vars)) {
    tvar  <- test.vars[i]
    tvec  <- data[subset, tvar, drop=FALSE]
    mflag <- any(is.na(tvec))
    if (mflag) {
      obj1 <- try(getAllObjects(data, cc.var, subtype.vars, tvar, covars.obj, options))
      tvec <- data[obj1$subset, tvar, drop=FALSE]
      tmp  <- try(fixedEffect.main(tvec, obj1, options))
    } else {
      tmp  <- try(fixedEffect.main(tvec, obj0, options))
    }
    ret[[tvar]] <- tmp
  }

  list(results=ret, stage1to2.model=obj0$z.all, outcomes=obj0$outcomes, 
       outcomes.removed=obj0$y.obj.rem)

} # END: fixedEffect

randomEffect <- function(data, cc.var, subtype.vars, test.vars, covars.obj=NULL, 
                 options=NULL) {

  options      <- check_options(options, fixed=FALSE)
  subtype.vars <- unique(subtype.vars)
  test.vars    <- unique(test.vars)
  checkDataAndVars(data, cc.var, subtype.vars, test.vars, covars.obj)

  ret <- list()
  for (i in 1:length(test.vars)) {
    tvar <- test.vars[i]
    options$snpv <- tvar
    obj1 <- try(getAllObjects(data, cc.var, subtype.vars, tvar, covars.obj, options))
    tvec <- data[obj1$subset, tvar, drop=FALSE]
    tmp  <- try(randomEffect.main(tvec, obj1, options))
    ret[[tvar]] <- tmp
  }

  list(results=ret, stage1to2.model=obj1$z.all, outcomes=obj1$outcomes, 
       outcomes.removed=obj1$y.obj.rem)

} # END: randomEffect


scoreTest <- function(x, support.obj, x.cov, z.design, debug=0) {

  debug                 <- as.integer(debug)
  zc                    <- support.obj$zc
  nparm.intere          <- as.integer(ncol(z.design))
  efficient.info        <- matrix(0,nparm.intere,nparm.intere)
  efficient.info.vec    <- as.numeric(efficient.info)
  info.complete.vec     <- efficient.info.vec
  info.lost.vec         <- efficient.info.vec
  score                 <- as.numeric(rep(0,nparm.intere))
  zc.nr                 <- as.integer(nrow(zc))
  zc.nc                 <- as.integer(ncol(zc))
  z.intere.nr           <- as.integer(nrow(z.design))
  z.intere.nc           <- as.integer(ncol(z.design))
  M                     <- z.intere.nr
  N                     <- as.integer(nrow(x))

  temp <- .C("ScoreTest",
               as.numeric(x[,1]),
               as.numeric(z.design),
               as.numeric(support.obj$inv_info_vec),
               as.numeric(support.obj$YminusP),
               as.numeric(support.obj$W_obs),
               score = score,
               efficient.info.vec=efficient.info.vec,
               zc.nr,
               zc.nc,
               z.intere.nr,
               z.intere.nc,
               nparm.intere,
               M,
               N,
               debug,
               info.complete.vec=info.complete.vec,
               info.lost.vec=info.lost.vec,
               as.numeric(x.cov),
               as.numeric(zc))
    score.result     <- matrix(temp$score, nrow=1)
    efficient.info   <- matrix(temp$efficient.info.vec,nparm.intere,nparm.intere)
    info.complete    <- matrix(temp$info.complete.vec,nparm.intere,nparm.intere)
    info.lost        <- matrix(temp$info.lost.vec,nparm.intere,nparm.intere)

  return(list(score.result=score.result, efficient.info.result=efficient.info, 
         info.complete=info.complete, info.lost=info.lost))

} # END: scoreTest

getSupport <- function(a, op) {

  # a :  list of all objects
  # op:  list of options

  z         <- a$z.all 
  newsnp    <- NULL
  fixed     <- op$fixed.effects
  if (!fixed) newsnp <- a[["snpv", exact=TRUE]]
  tmp       <- getDelta0ForSupport(a, op, newsnp=newsnp)
  delta0    <- tmp$delta0
  nuis.rows <- a$nuis.rows.1
  nuis.cols <- a$nuis.cols.2
  zc        <- z[nuis.rows, nuis.cols, drop=FALSE]

  if ((op$fixed.effects) || (is.null(newsnp))) {
    z         <- zc
    delta0    <- delta0[nuis.cols]
    x         <- a$x.obj
  } else {
    # Apply the subset here
    x         <- cbind(newsnp[a$subset], a$x.obj)
  }

  ret <- scoreSupport(a$prob.fit.result, x, z, delta0, a$M, tol=op$tol, 
                      tolMaxstep=op$tolMaxstep, NITER=op$maxiter, debug=op$debug) 
  ret$add.snp  <- tmp$add.snp
  ret$flip.snp <- tmp$flip.snp
  
  # For score test, zc should not contain the test rows/cols
  ret$zc <- zc

  ret

} # END: getSupport

scoreSupport <- function(prob.fit.result, x.all, z.all, delta0, M, 
                               tol=1e-4, tolMaxstep=1e-3, NITER=100, debug=0) {

  tol              <- as.numeric(tol)
  tolMaxstep       <- as.numeric(tolMaxstep)
  N                <- as.integer(nrow(x.all))
  NCOV             <- as.integer(ncol(x.all))
  NM               <- N*M
  nparm            <- as.integer(length(delta0))
  NITER            <- as.integer(NITER)
  Znr              <- as.integer(nrow(z.all))
  Znc              <- as.integer(ncol(z.all))
  debug            <- as.integer(debug)
  ret_rc           <- as.integer(1)
  ret_delta        <- as.numeric(rep(-9999, nparm))
  ret_info         <- as.numeric(rep(-9999,nparm^2))
  ret_p            <- as.numeric(rep(0,NM))
  ret_Inv_info_vec <- as.numeric(as.vector(matrix(0,Znc,Znc)))
  YminusP          <- as.numeric(prob.fit.result$y_em)
  W_obs            <- as.numeric(rep(0,NM*M))
  missing.number   <- as.integer(prob.fit.result$missing.number)

  temp <- .C("ScoreSupport",
             as.numeric(delta0),
             nparm,
             as.numeric(prob.fit.result$y_em),
             as.numeric(x.all),
             as.numeric(z.all),
             Znr,Znc, N, M, NCOV, NITER,
             tol,
             tolMaxstep,
             debug,
             ret_rc=ret_rc,
             ret_delta=ret_delta,
             ret_info=ret_info,
             ret_p=ret_p,
             as.integer(prob.fit.result$missing.vec),
             as.integer(prob.fit.result$missing.mat),
             missing.number,
             ret_Inv_info_vec=ret_Inv_info_vec,
             YminusP=YminusP,
             W_obs = W_obs)
  #print(paste0("EM Algorithm Converged"))
  inv_info_vec <- temp$ret_Inv_info_vec
  YminusP      <- temp$YminusP
  W_obs        <- temp$W_obs
  ret_delta    <- temp$ret_delta

#print(paste("max_abs_diff=", max(abs(ret_delta - delta0), sep="")))
#print(abs(delta0-ret_delta))
#print("######")

  return(list(inv_info_vec=inv_info_vec,YminusP=YminusP,W_obs=W_obs, delta=ret_delta))
  
} # END: scoreSupport

fixedEffect.main <- function(x, obj, op, add.parm.names=TRUE) {

  score.obj  <- scoreTest(x, obj$support.obj, obj$x.obj, obj$z.design, debug=op$debug)
  score.test <- ScoreGlobalTestForAssoc2(score.obj)

  ret <- list(global.score.test=score.test, score.estimates=score.obj[["score.result", exact=TRUE]],
              score.information=score.obj[["efficient.info.result", exact=TRUE]],
              stage2.estimates=obj$support.obj$delta)

  # Add parameter names to objects
  if (add.parm.names) {
    tmp <- c("score.estimates", "score.information", "stage2.estimates") 
    ret <- addParmNames(obj, ret, tmp, op)
  }

  ret

} # END: fixedEffect.main

randomEffect.main <- function(x, obj, op, add.parm.names=TRUE) {

  score.obj  <- scoreTest(x, obj$support.obj, obj$x.obj, obj$z.design, debug=op$debug)
  score.test <- ScoreGlobalMixedTestForAssoc2(score.obj)

  ret <- list(global.score.test=score.test, score.estimates=score.obj[["score.result", exact=TRUE]],
              score.information=score.obj[["efficient.info.result", exact=TRUE]],
              stage2.estimates=obj$support.obj$delta)

  # Add parameter names to objects
  if (add.parm.names) {
    tmp <- c("score.estimates", "score.information", "stage2.estimates") 
    ret <- addParmNames(obj, ret, tmp, op)
  }

  ret

} # END: randomEffect.main

mixedEffect.main <- function(which, x, obj, op, add.parm.names=TRUE) {

  if (is.null(ncol(x))) dim(x) <- c(1, length(x))
  if (which == 1) {
    ret <- fixedEffect.main(x, obj, op, add.parm.names=add.parm.names)
  } else {
    ret <- randomEffect.main(x, obj, op, add.parm.names=add.parm.names)
  }

  ret

} # END: mixedEffect.main

# Function to get model info
modelInfo <- function(data, cc.var, subtype.vars, covars.obj=NULL, 
                 options=NULL) {

  options      <- check_options(options)
  subtype.vars <- unique(subtype.vars)
  checkDataAndVars(data, cc.var, subtype.vars, NULL, covars.obj, check.test=FALSE)

  obj <- getAllObjects(data, cc.var, subtype.vars, NULL, covars.obj, options, retPos=1)

  valid <- getValidOptions()$valid
  all   <- names(options)
  rem   <- all[!(all %in% valid)]
  if (length(rem)) {
    for (tmp in rem) options[[tmp]] <- NULL
  } 
  options$test.cols <- obj$test.cols.2

  list(stage1to2.model=obj$z.all, outcomes=obj$outcomes, 
       outcomes.removed=obj[["y.obj.rem", exact=TRUE]],
       delta0=obj$delta0, options=options)

} # END: getModelInfo
