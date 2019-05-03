
#' Title
#'
#' @param delta0 
#' @param y 
#' @param x.all 
#' @param z.standard 
#' @param z.all 
#'
#' @return
#' @export
#'
#' @examples
StepScoreTestSupportMixedModel <- function(delta0,y,x.all,z.standard,z.all){
  tol <- as.numeric(1e-04)
  delta_old <- delta0
  prob.fit.result <- ProbFitting(delta0,y,x.all,z.standard,z.all,missingTumorIndicator=NULL)
    y.fit <- prob.fit.result[[1]]
    idx.drop = prob.fit.result[[4]]
    if(length(idx.drop)!=0){
      x.all <- x.all[-idx.drop,,drop=F]
      y_em <- y_em[-idx.drop,,drop=F]
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
    Y <- as.numeric(as.vector(y.fit))
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
    
    
    
    
    temp <- .C("CompleteCasesScoreTestSupport",
               deltai,
               nparm,
               Y=Y,
               X,
               ZallVec,
               Znr,
               Znc,
               N,
               M,
               NCOV,
               NITER,
               tol,
               debug,
               ret_rc=ret_rc,
               ret_delta=ret_delta,
               ret_info=ret_info,
               ret_p=ret_p,
               ret_Inv_info_vec=ret_Inv_info_vec,
               YminusP=YminusP,
               W= W_obs,
               WXZ_vec = WXZ_vec,
               WX_vec = WX_vec)
    print(paste0("Algorithm Converged"))
    inv_info_vec <- temp$ret_Inv_info_vec
    YminusP <- temp$YminusP
    W_obs <- temp$W
    WXZ_vec <- temp$WXZ_vec
    WX_vec <- temp$WX_vec
    
    
    result <- list(inv_info_vec=inv_info_vec,YminusP=YminusP,W_obs=W_obs,WXZ_vec = WXZ_vec,zc=z.all)
    return(result)
  }