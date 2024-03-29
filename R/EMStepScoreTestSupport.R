#' Title
#'
#' @param delta0
#' @param y
#' @param x.all
#' @param z.standard
#' @param z.all
#' @param missingTumorIndicator
#'
#' @return
#' @export
#'
#' @examples
EMStepScoreTestSupport <- function(delta0,y,x.all,z.standard,z.all,missingTumorIndicator){

  tol <- as.numeric(1e-04)
  tolMaxstep <- as.numeric(1e-03)
  #delta_old <- rep(0,length(delta0))
  delta_old <- delta0
  ##EM algorithm
  ##first E step
  #print(paste0("Begin EM algorithm"))
  #print(paste0("EM round: 1"))
  prob.fit.result <- ProbFitting(delta_old,as.matrix(y),x.all,
                                 z.standard,z.all,missingTumorIndicator)
  y_em <- prob.fit.result[[1]]
  missing.vec <- as.numeric(as.vector(prob.fit.result[[2]]))
  missing.mat <- prob.fit.result[[3]]
  missing.mat.vec <- as.numeric(as.vector(missing.mat))
  missing.number <- as.integer(length(missing.vec))
  idx.drop = prob.fit.result[[4]]
  #remove all the subtypes below the cutoff threshold
  if(length(idx.drop)!=0){
    x.all <- x.all[-idx.drop,,drop=F]
    y_em <- y_em[-idx.drop,,drop=F]
  }
  for(k in 1:length(missing.vec)){
    missing.vec[k] <- missing.vec[k]-sum(missing.vec[k]>=idx.drop)
  }
  
  # sof <- "try5.so"
  # dyn.load(sof)

  N <- as.integer(nrow(x.all))
  #x <- cbind(1,x)
  #p <- ncol(x)
  M <- as.integer(nrow(z.standard))

  NCOV   <- as.integer(ncol(x.all))
  NM     <- N*M
  nparm  <- as.integer(length(delta0))
  deltai <- as.numeric(delta0)

  NITER  <- as.integer(500)
  Y <- as.numeric(as.vector(y_em))
  X <- as.numeric(as.vector(x.all))
  ZallVec = as.numeric(as.vector(z.all))
  Znr = as.integer(nrow(z.all))
  Znc = as.integer(ncol(z.all))
  debug     <- as.integer(1)
  ret_rc    <- as.integer(1)
  ret_delta <- as.numeric(rep(-9999, nparm))
  ret_info <- as.numeric(rep(-9999,nparm^2))
  ret_p <- as.numeric(rep(0,NM))
  ret_Inv_info_vec <- as.numeric(as.vector(matrix(0,Znc,Znc)))
  YminusP <- Y
  W_obs <- as.numeric(rep(0,N*M*M))
  WXZ_vec <- as.numeric(rep(0,N*M*Znc))
  WX_vec <- as.numeric(rep(0,N*M*Znr))


  temp <- .C("EMStepScoreSupport",
             deltai,
             nparm,
             Y=Y,
             X,
             ZallVec,
             Znr,Znc, N, M, NCOV, NITER,
             tol,
             tolMaxstep,
             debug,
             ret_rc=ret_rc,
             ret_delta=ret_delta,
             ret_info=ret_info,
             ret_p=ret_p,
             missing.vec,
             missing.mat.vec,missing.number,
             ret_Inv_info_vec=ret_Inv_info_vec,
             YminusP=YminusP,
             W_obs = W_obs,
             WXZ_vec = WXZ_vec,
             WX_vec = WX_vec)
  print(paste0("EM Algorithm Converged"))
  inv_info_vec <- temp$ret_Inv_info_vec
  YminusP <- temp$YminusP
  W_obs <- temp$W_obs
  WXZ_vec <- temp$WXZ_vec
  WX_vec <- temp$WX_vec


  return(list(inv_info_vec=inv_info_vec,YminusP=YminusP,W_obs=W_obs,WXZ_vec = WXZ_vec,zc=z.all,idx.drop = idx.drop))
  #  return(temp)
}

