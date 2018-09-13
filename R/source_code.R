
check_delta0 <- function(d, nc=0) {

  if (!is.numeric(d)) stop("ERROR: options$delta0 must be an numeric vector")
  if ((nc) && (length(d) != nc)) {
    stop("ERROR: length(options$delta0) must be equal to ncol(options$stage1to2.model)")
  }

} # END: check_delta0

check_model <- function(op) {

  v   <- c("additive", "baselineonly", "pairwise.interaction", "saturated")
  mod <- op$stage1to2.model
  if (is.character(mod) && (length(mod) == 1)) {
    if (!(mod %in% v)) {
      cat("\n ERROR options$model must be a matrix or one of the following:\n")
      cat(v)
      stop("ERROR: with options$model")
    }
    op$test.cols <- NULL
  } else if (is.matrix(mod)) {
    # In this case test.cols must be specified
    nc <- ncol(mod)
    tc <- op[["test.cols", exact=TRUE]]
    if (is.null(tc)) stop("ERROR: options$test.cols must be specified when options$stage1to2 is a matrix")
    if (!is.numeric(tc)) stop("ERROR: options$test.cols must be an integer vector")
    if (any(tc > nc)) stop("ERROR: with options$test.cols")

    # Check delta0
    d <- op[["delta0", exact=TRUE]]
    if (!is.null(d)) check_delta0(d, nc=nc) 
  } else {
    stop("ERROR: with options$model")
  }

  op

} # END: check_model

getValidOptions <- function() {

  valid  <- c("stage1to2.model", "tol", "tolMaxstep", "debug", "maxiter",
               "min.ncase", "min.ncontrol", "delta0", "test.cols") 
  values <- list("additive", 1e-4, 1e-4, 0, 100,
                  10, 10, NULL, NULL)

  list(valid=valid, values=values)

} # END: getValidOptions

check_options <- function(op, fixed=TRUE) {

  tmp     <- getValidOptions()
  valid   <- tmp$valid
  values  <- tmp$values
  opNames <- names(op)
  tmp     <- !(opNames %in% valid)
  m       <- sum(tmp)
  if (m) {
    str <- paste(opNames[tmp], collapse=", ", sep="")
    if (m > 1) {
      str <- paste("ERROR: ", str, " are not valid options", sep="")
    } else {
      str <- paste("ERROR: ", str, " is not a valid option", sep="")
    }
    stop(str)
  }
  op  <- default.list(op, valid, values)
  op  <- check_model(op) 
  
  if (op$tol <= 0) stop("ERROR: options$tol must be positive")
  if (op$tolMaxstep <= 0) stop("ERROR: options$tolMaxstep must be positive")
  if (op$min.ncase <= 0) stop("ERROR: options$min.ncase must be positive")
  if (op$min.ncontrol <= 0) stop("ERROR: options$min.ncontrol must be positive")
  if (op$maxiter <= 0) stop("ERROR: options$maxiter must be positive")
  
  op$fixed.effects      <- fixed
  op$miss.ind           <- 888
  op$stopOnDelta0.error <- TRUE

  op

} # END: check_options

is.formula <- function(x){
   inherits(x,"formula")
} # END: is.formula

getVars <- function(obj) {

  if (is.null(obj)) return(NULL)
  if (is.list(obj)) {
    vars <- unique(unlist(obj))
  } else if (is.formula(obj)) {
    vars <- unique(all.vars(obj)) 
  } else {
    if (is.vector(obj)) {
      vars <- obj
    } else {
      stop("ERROR with obj")
    }
  }

  vars

} # END: getVars

checkModelVars <- function(cc.var, subtype.vars, test.vars, covars.obj, check.test=TRUE) {

  subFlag <- !is.null(subtype.vars)
  if (length(subtype.vars) < 2) stop("ERROR: There must be at least 2 subtype.vars") 
  if (length(cc.var) != 1) stop("ERROR: cc.var must be a single variable name") 
  if (check.test) {
    if (length(test.vars) < 1) stop("ERROR: test.vars is empty") 
    if (cc.var %in% test.vars) stop("ERROR: cc.var must not be in test.vars")
  }
  covFlag <- !is.null(covars.obj)
  covars  <- getVars(covars.obj)
  if (subFlag) {
    if (length(subtype.vars) < 2) stop("ERROR: currently, there must be at least 2 subtype.vars") 
    if (cc.var %in% subtype.vars) stop("ERROR: subtype.vars must not contain cc.var")
  }
  if (covFlag) {
    if ((!is.formula(covars.obj)) && (!is.vector(covars.obj))) stop("ERROR: covars must be a vector or formula") 
    if (cc.var %in% covars) stop("ERROR: covars must not contain cc.var")
    if ((check.test) && (test.vars %in% covars)) stop("ERROR: covars must not contain test.vars")
    if (subFlag) {
      if (any(subtype.vars %in% covars)) stop("ERROR: covars must not contain any subtype.vars")
    }
  }

  covars

} # END: checkModelVars

checkDataVars <- function(data, vars, bin.vars=NULL, num.vars=NULL) {

  vars <- unique(c(vars, bin.vars, num.vars))
  tmp  <- !(vars %in% colnames(data))
  if (any(tmp)) {
    miss <- vars[tmp]
    tmp  <- paste(miss, collapse=", ", sep="")
    print("The following variables were not found in data:")
    print(tmp)
    stop("ERROR in checkDataVars")
  }

  # Binary variables
  if (!is.null(bin.vars)) {
    for (var in bin.vars) {
      vec <- unique(data[, var])
      if (!all(vec %in% c(0, 1, NA))) {
        str <- paste("ERROR: the variable ", var, " is not binary", sep="")
        stop(str)
      }
    }
  }

  # Numeric variables
  if (!is.null(num.vars)) {
    for (var in num.vars) {
      if (!is.numeric(data[, var])) {
        str <- paste("ERROR: the variable ", var, " is not numeric", sep="")
        stop(str)
      }
    }
  }

  NULL

} # END: checkDataVars

checkDataAndVars <- function(data, cc.var, subtype.vars, test.vars, covars.obj,
                             check.test=TRUE) {
  
  if (!is.data.frame(data)) stop("ERROR: data must be a data frame")
  covars   <- checkModelVars(cc.var, subtype.vars, test.vars, covars.obj, check.test=check.test) 
  allvars  <- unique(c(cc.var, subtype.vars, test.vars, covars))
  num.vars <- c(cc.var, subtype.vars, test.vars)
  if (!is.formula(covars.obj)) num.vars <- c(num.vars, covars)
  binvars <- NULL
  if (!is.null(subtype.vars)) binvars=cc.var
  checkDataVars(data, allvars, bin.vars=binvars, num.vars=num.vars)

  NULL
 
} # END: checkDataAndVars

getYobj <- function(data, cc.var, subtype.vars, subset, op) {

  yvars <- c(cc.var, subtype.vars)
  y.obj <- as.matrix(data[subset, yvars, drop=FALSE])

  tmp   <- !is.finite(y.obj)

  if (any(tmp)) y.obj[tmp] <- op$miss.ind
  
  # Subtype vars should be coded as NA for controls
  tmp <- y.obj[, 1] %in% 0
  y.obj[tmp, 2:ncol(y.obj)] <- NA

  # Check the number of controls
  if (sum(tmp) < op$min.ncontrol) stop("ERROR: too few controls in data") 

  # Check the number of cases
  ret <- checkYobjFreq(y.obj, op)

  ret

} # END: getYobj

getCompleteCases <- function(y.obj) {

  mat <- is.na(y.obj)
  vec <- rowSums(mat)
  tmp <- vec > 0
  return(!tmp)

} # END: getCompleteCases

getYobjFreq <- function(y.obj) {

  nc <- ncol(y.obj)
  if (nc < 3) stop("ERROR with y.obj")

  # Get strings to identify subtypes
  vec <- y.obj[, 2]
  for (i in 3:nc) vec <- paste(vec, y.obj[, i], sep=".")

  # Get the complete cases
  case <- getCompleteCases(y.obj)
  
  # Get the frequency counts among the cases
  freqs <- table(vec[case])

  list(vec=vec, freqs=freqs)

} # END: getYobjFreq

checkYobjFreq <- function(y.obj, op) {

  tmp     <- getYobjFreq(y.obj)
  freqs   <- tmp$freqs 
  vec     <- tmp$vec  
  removed <- NULL

  # Remove subtypes with too few cases
  tmp   <- freqs < op$min.ncase
  if (any(tmp)) {
    nn      <- names(freqs)
    rem     <- nn[tmp]
    tmp     <- vec %in% rem
    removed <- unique(y.obj[tmp, -1, drop=FALSE])
    rownames(removed) <- NULL
    y.obj   <- y.obj[!tmp, , drop=FALSE]
  }

  # Check cases
  case <- y.obj[, 1] %in% 1
  if (!sum(case)) stop("ERROR: all cases have been removed")

  list(y.obj=y.obj, removed=removed)

} # END: checkYobjFreq

getInitialSubset <- function(data, cc.var, test.vars) {

  subset <- is.finite(data[, cc.var])

  # Remove the common missing values from all test variables
  if (length(test.vars)) {
    vec    <- rowSums(!is.finite(as.matrix(data[, test.vars, drop=FALSE])))
    tmp    <- vec == length(test.vars)
    subset <- subset & !tmp
    subset[is.na(subset)] <- FALSE
  }

  if (!sum(subset)) stop("ERROR: all rows have been removed")

  subset

} # END: getInitialSubset

getXobj <- function(data, covars.obj, subset, op) {

  # The returned object will be the design matrix of covariates
  mat <- NULL
  if (is.null(covars.obj)) return(NULL)

  if (is.formula(covars.obj)) {
    mat    <- model.matrix(covars.obj, data=data)
  } else {
    covars <- getVars(covars.obj)
    mat    <- as.matrix(data[, covars, drop=FALSE]) 
    rownames(mat) <- rownames(data)
  }
  if (!nrow(mat)) stop("ERROR with covars.obj: all rows in data have been removed")

  vec <- rowSums(!is.finite(mat))
  tmp <- vec > 0
  if (any(tmp)) mat <- mat[!tmp, , drop=FALSE]
  if (!nrow(mat)) stop("ERROR with covars.obj: all rows in data have been removed")

  mat  

} # END: getXobj

getDataObjects <- function(data, cc.var, subtype.vars, test.vars, covars.obj, op) {

  nr  <- nrow(data)
  ids <- rownames(data)
  if (is.null(ids)) {
    # Should never be NULL to begin with
    ids <- paste("r", 1:nr, sep="")
    rownames(data) <- ids
  }
  covFlag <- !is.null(covars.obj)
  mat <- NULL

  # Get the subset of subs to use initially
  subset <- getInitialSubset(data, cc.var, test.vars)

  # Get the design matrix of covariates  
  if (covFlag) {
    mat <- getXobj(data, covars.obj, subset, op)

    # Update subset
    tmp <- !(rownames(data) %in% rownames(mat))
    if (any(tmp)) subset[tmp] <- FALSE
  }
  
  # Get the y matrix
  tmp       <- getYobj(data, cc.var, subtype.vars, subset, op)
  y.obj     <- tmp$y.obj 
  y.obj.rem <- tmp$removed
 
  # Update subset and mat
  tmp <- !(rownames(data) %in% rownames(y.obj))
  if (any(tmp)) subset[tmp] <- FALSE
  if (covFlag) mat <- mat[subset, , drop=FALSE]

  # Remove any constant columns (including intercept)
  if (covFlag) {
    tmp <- rep(TRUE, ncol(mat))
    for (i in 1:ncol(mat)) {
      if (length(unique(mat[, i])) < 2) tmp[i] <- FALSE
    }
    mat <- mat[, tmp, drop=FALSE]
    if (!ncol(mat)) {
      warning("All covariates have been removed")
      mat <- NULL
    }
  }

  list(y.obj=y.obj, x.obj=mat, subset=subset, y.obj.rem=y.obj.rem)

} # END: getDataObjects

setZdesignMat <- function(obj, op) {

  model   <- op$stage1to2.model
  userDef <- obj$userDefined
  if (!userDef) {
    str   <- paste("z.design.", model, sep="")
    mat   <- as.matrix(obj[[str, exact=TRUE]])
  } else {
    rows  <- obj$test.rows.1
    cols  <- obj$test.cols.2 
    mat   <- obj$z.all[rows, cols]
  }
  
  obj$z.design.additive              <- NULL
  obj$z.design.baselineonly          <- NULL
  obj$z.design.pairwise.interaction  <- NULL
  obj$z.design.saturated             <- NULL
  obj$z.deisign.baselineonly         <- NULL
  obj$z.design                       <- mat
  obj$z.standard                     <- NULL

  obj

} # END: setZdesignMat

# For the row names of z.all (not user-defined)
getParmNames.stage1 <- function(outcomes, varnames) {

  # varnames should not include intercept or the test variable
  outcomes <- sort(outcomes)
  ret      <- NULL
  for (v in outcomes) ret <- c(ret, paste(varnames, v, sep="."))
  
  ret

} # END: getParmNames.stage1

# For column names of z.all and delta (not user-defined)
getParmNames.stage2 <- function(n, test.cols=NULL) {

  ret <- paste("Delta", 1:n, sep="")
  if (!is.null(test.cols)) ret[test.cols] <- paste(ret[test.cols], "(Test)", sep="")
  
  ret

} # END: getParmNames.stage2


getOutcomeMap <- function(a) {

  M                  <- a$M
  outcomes           <- a$freq.subtypes
  nc                 <- ncol(outcomes)
  tmp                <- outcomes[, nc] > 0
  tmp[is.na(tmp)]    <- FALSE
  outcomes           <- outcomes[tmp, -nc, drop=FALSE]
  n                  <- nrow(outcomes)
  outcomes           <- cbind(1:n, outcomes)
  ynames             <- colnames(a$y.obj)
  colnames(outcomes) <- c("Outcome", ynames[-1]) 

  outcomes

} # END: getOutcomeMap

getDefaultDelta0 <- function(a) {

  # a is the list of all objects
  
   z     <- a$z.all
   delta <- rep(0, ncol(z))
   names(delta) <- colnames(z)
   
   if (!(a$userDefined)) { 
     # Deltas for intercepts have been computed already
     M            <- a$M
     d0           <- a$delta0
     delta[1:M]   <- d0[1:M]
   }

   delta

} # END: getDefaultDelta0

getDelta0 <- function(a, op) {

  np <- ncol(a$z.all)

  # If a user specified delta0, then use it provided it has the correct length
  delta <- op[["delta0", exact=TRUE]]
  if (!is.null(delta)) {
    tmp <- try(check_delta0(delta, nc=np), silent=TRUE)
    if (checkTryError(tmp)) {
      if (op$stopOnDelta0.error) {
        stop(tmp)
      } else {
        delta <- getDefaultDelta0(a)
      }
    }
  } else {
    delta <- getDefaultDelta0(a)
  }

  delta

} # END: getDelta0

getDelta.R2 <- function(newsnp, deltaList) {

    deltaN   <- deltaList$deltaN
    corr     <- cor(deltaList$snpMat[, 1:deltaN, drop=FALSE], newsnp, method="pearson", use="pairwise.complete.obs")
    corr     <- as.vector(corr) 
    tmp      <- !is.finite(corr)
    if (any(tmp)) corr[tmp] <- 0
    r2       <- corr*corr
    jj   <- which.max(r2)
    ret  <- deltaList$deltaMat[, jj]
    sgn  <- sign(corr[jj])
    flip <- sgn < 0
    if (flip) {
      # If the correlation is negative, flip signs of the snp parms
      ids      <- deltaList$test.cols.2
      ret[ids] <- -ret[ids]
    }
    add <- (all(r2 < deltaList$max.r2)) 

    list(delta=ret, add.snp.flag=add, flip.snp=flip, maxR2=r2[jj])
  
} # END: getDelta.R2

# Function to get initial estimates for support function
getDelta0ForSupport <- function(a, op, newsnp=NULL) {

  if (is.null(op[["snpv", exact=TRUE]])) {
    scan    <- FALSE
  } else {
    scan    <- TRUE
  }

  # Get the length of delta0 (all parms included)
  np       <- ncol(a$z.all)
  cols2    <- a$nuis.cols.2
  add.snp  <- FALSE
  flip.snp <- FALSE

  # First check for a user defined delta0
  delta0 <- op[["delta0", exact=TRUE]]
  
  if (is.null(delta0)) {
    # Not user-defined

    fixed <- op$fixed.effects
    if (fixed) {
      # Use the fitted parms (if they exist)
      del.nuis <- op[["delta.nuis", exact=TRUE]]

      if (length(del.nuis) == length(cols2)) {
        delta0        <- rep(0, np)
        delta0[cols2] <- del.nuis
      } else {
        delta0        <- getDefaultDelta0(a)
      }
    } else {
      # For random effects, try using R2 first
      snpFlag <- !is.null(newsnp)
      dlist   <- op[["delta0.list", exact=TRUE]]
      if ((!is.null(dlist)) && (snpFlag)) {
        if (dlist$deltaN) {
          tmp      <- getDelta.R2(newsnp, dlist)
          delta0   <- tmp$delta
          add.snp  <- tmp$add.snp
          flip.snp <- tmp$flip.snp
        } else {
          # An initial snp must be added
          add.snp  <- TRUE
        }
      }
      if (is.null(delta0)) {
        # Use the fitted parms (if they exist)
        del.nuis <- op[["delta.nuis", exact=TRUE]]
        if (length(del.nuis) == length(cols2)) {
          delta0        <- rep(0, np)
          delta0[cols2] <- del.nuis
        } else {
          delta0        <- getDefaultDelta0(a)
        }
      }
    }
  }

  # Final check of delta0
  tmp <- try(check_delta0(delta0, nc=np), silent=TRUE)
  if (checkTryError(tmp)) {
    if (!scan) {
      stop(tmp)
    } else {
      delta0 <- getDefaultDelta0(a)
    }
  }

  list(delta0=delta0, add.snp=add.snp, flip.snp=flip.snp)

} # END: getDelta0ForSupport

# Function to return z.all matrix
getZ <- function(a, op) {

  model   <- op$stage1to2.model
  userDef <- a$userDefined
  M       <- a$M
  ncov    <- a$ncov
  ncovp1  <- ncov + 1
  ncovp2  <- ncov + 2
  NR      <- M*ncovp2 # Add 1 for intercept and test parms
  
  # Get row and column numbers for test var and nuisance vars
  int.rows <- seq(from=1, to=NR, by=ncovp2)
  int.cols <- 1:M

  # Stage 1 rows in z.all for test variable (variable of interest)
  test.rows.1 <- int.rows + 1
  nuis.rows.1 <- (1:NR)[-test.rows.1]

  if (!userDef) {
    str <- paste("z.design.", model, sep="")
    mat <- a[[str, exact=TRUE]]
    znc <- ncol(mat)
    NC  <- M + ncovp1*znc  # Add 1 for test var

    # Stage 2 cols in z.all for test variable (variable of interest)
    test.cols.2 <- (M+1):(M+znc)
  } else {
    z.all <- model
    # Check the number of rows
    if (nrow(z.all) != NR) stop("ERROR: incorrect number of rows for options$stage1to2.model")
    NC          <- ncol(z.all)
    test.cols.2 <- op$test.cols # This must be specified by the user
  }
  
  # Indices for nuisance parms
  nuis.cols.2 <- (1:NC)[-test.cols.2]

  if (!userDef) {
    z.all <- matrix(data=0, nrow=NR, ncol=NC)
    for (i in 1:M) z.all[int.rows[i], int.cols[i]] <- 1
    for (j in 1:M) {
      vec <- as.numeric(mat[j, ])
      bb  <- M
      for (i in 1:ncovp1) {
        row   <- (j-1)*ncovp2 + i + 1
        aa    <- bb + 1
        bb    <- aa + znc - 1
        cols  <- aa:bb 
        z.all[row, cols] <- vec
      }
    }

    cx              <- NULL
    if (ncov) cx    <- colnames(a$x.obj)
    tmp             <- c("(Intercept)", "(Test)", cx)
    parms.stage1    <- getParmNames.stage1(1:M, tmp)
    parms.stage2    <- getParmNames.stage2(NC)
    rownames(z.all) <- parms.stage1
    colnames(z.all) <- parms.stage2
  }

  list(z.all=z.all, test.rows.1=test.rows.1, nuis.rows.1=nuis.rows.1,
       test.cols.2=test.cols.2, nuis.cols.2=nuis.cols.2)

} # END: getZ

# Function to update objects depending on a fixed or random effects analysis
updateObjects <- function(a, op) {

    test.rows <- a$test.rows.1
    test.cols <- a$test.cols.2
    nuis.rows <- a$nuis.rows.1
    nuis.cols <- a$nuis.cols.2
    z.all     <- a$z.all
    delta0    <- a$delta0
    a$z.all   <- z.all[nuis.rows, nuis.cols, drop=FALSE]
    a$delta0  <- delta0[nuis.cols]

    a

} # END: updateObjects

getMissingMat <- function(y.obj, z.standard, missingTumorIndicator) {

  y1         <- y.obj[, 1]
  y          <- y.obj[, -1, drop=FALSE]
  nc         <- ncol(y)
  nr         <- nrow(y)
  tmpy1      <- y1 == 1
  tmp        <- y == missingTumorIndicator
  notMissMat <- !tmp
  sumMiss    <- rowSums(tmp)
  miss       <- (sumMiss > 0) & tmpy1
  miss[is.na(miss)] <- FALSE
  nmiss <- sum(miss)
  if (!nmiss) return(list(missing.vec=NULL, missing.mat=NULL)) 

  missing.vec <- (1:nr)[miss] 
  M           <- nrow(z.standard)
  missing.mat <- matrix(0,nrow=nmiss, ncol=M)
  oneToM      <- 1:M
  for (i in 1:nmiss) {
    row  <- missing.vec[i]
    idx  <- notMissMat[row, ]
    nidx <- sum(idx)
    if (!nidx) {
      missing.mat[i, ] <- 1
    } else {
      # Find all rows of z.standard that match
      tmp <- z.standard[, idx, drop=FALSE] == matrix(y[row, idx], byrow=TRUE, nrow=M, ncol=nidx)  
      tmp <- rowSums(tmp) == nidx
      if (any(tmp)) {
        missing.mat[i, ] <- as.numeric(tmp)
      } else {
        missing.mat[i, ] <- 1 
      }
    }
  }

  list(missing.vec=missing.vec, missing.mat=missing.mat)

} # END: getMissingMat

# Function to get all objects
getAllObjects <- function(data, cc.var, subtype.vars, test.vars, covars.obj, op, retPos=0) {

  # Get the Y, X and test objects
  obj <- getDataObjects(data, cc.var, subtype.vars, test.vars, covars.obj, op)

  obj$ncov <- CountCovarNumber(obj[["x.obj", exact=TRUE]])
  obj$userDefined <- length(op$stage1to2.model) > 1

  y.obj     <- obj$y.obj
  case      <- !(y.obj[, 1] %in% 0)
  missInd   <- op$miss.ind
  missFlag  <- FALSE
  fixedFlag <- op$fixed.effects

  # Missing values have already been changed to op$missInd
  if (any(y.obj[case, ] %in% missInd)) missFlag <- TRUE
  obj$missFlag <- missFlag

  if (missFlag) {
    missing.data.vec <- GenerateMissingPosition(y.obj, missInd)
    y.pheno.complete <- y.obj[-missing.data.vec, , drop=FALSE]
  } else {
    y.pheno.complete <- y.obj
    missInd          <- NULL
  }  
  obj$missInd <- missInd

  # Add a test column to design matrix
  tmp   <- cbind(rep(1, nrow(y.obj)), obj$x.obj)
  cx    <- c("(Test)", colnames(obj$x.obj))
  colnames(tmp) <- cx

  # InitialSetup returns: list(delta0 = delta0,z.all=z.all,z.standard= z.standard,
  #  z.design.baselineonly = z.design.baselineonly,z.design.additive=z.design.additive,
  #  z.design.pairwise.interaction=z.design.pairwise.interaction,
  #  z.design.saturated=z.design.saturated, covar.names = covar.names,tumor.names=tumor.names))
  tmp <- InitialSetup(y.pheno.complete, NULL, tmp, NULL,
                      NULL, x.self.design=NULL, z.design=NULL, cutoff=op$min.ncase-1)
  obj   <- c(obj, tmp)
  obj$M <- nrow(obj$z.standard)
  M     <- obj$M

  # Get the full Z matrix that maps stage 1 to stage 2
  obj$z.all  <- NULL
  tmp        <- getZ(obj, op)
  obj        <- c(obj, tmp)
  obj$delta0 <- getDelta0(obj, op)

  # Check for errors
  ndelta <- length(obj$delta0)
  n1     <- nrow(obj$z.all)
  if (ndelta != ncol(obj$z.all)) stop("ERROR with z.all and/or delta0")
  if (ndelta > n1) stop("ERROR: more second stage parameters than first")
  if (ndelta == n1) stop("ERROR: same number of second stage parameters as first")

  # Get outcomes
  obj$outcomes <- getOutcomeMap(obj) 

  # If only certain info is needed return early
  if (retPos) return(obj)

  # y_em, missing.vec and missing.mat are needed. complete.vec is not needed
  tmp <- ProbFitting(obj$delta0[obj$nuis.cols.2], as.matrix(obj$y.obj), obj$x.obj, 
          obj$z.standard, obj$z.all[obj$nuis.rows.1, obj$nuis.cols.2], missInd)
  tmp$complete.vec <- NULL
  missing.vec      <- tmp[["missing.vec", exact=TRUE]]
  missing.number   <- length(missing.vec)
  if (!missing.number) {
    # Set these to something to pass into the C code
    tmp$missing.vec <- 0
    tmp$missing.mat <- 0
  } 
  tmp$missing.number <- missing.number
  obj$prob.fit.result <- tmp

  obj <- setZdesignMat(obj, op)

  # Save the snp if a random effects scan is running
  # We need all genotypes (including missing) for R^2
  # Apply the subset later in getScore
  if (!fixedFlag) {
    snpv <- op[["snpv", exact=TRUE]]
    if (!is.null(snpv)) obj$snpv <- data[, snpv]
  }

  # Get the score support
  obj$support.obj <- getSupport(obj, op)

  # Save original fitted (nuisance) parms
  obj$delta.nuis <- obj$support.obj$delta

  obj

} # END: getAllObjects

# Function to assign a default value to an element in a list
default.list <- function(inList, names, default, error=NULL,
                         checkList=NULL) {

  # inList      List
  # names       Vector of names of items in inList
  # default     List of default values to assign if a name is not found
  #             The order of default must be the same as in names.
  # error       Vector of TRUE/FALSE if it is an error not to have the
  #             name in the list. 
  #             The default is NULL
  # checkList   List of valid values for each name.
  #             Use NA to skip a list element.
  #             The default is NULL

  n1 <- length(names)
  n2 <- length(default)
  if (n1 != n2) stop("ERROR: in calling default.list")

  if (is.null(error)) {
    error <- rep(0, times=n1)
  } else if (n1 != length(error)) {
    stop("ERROR: in calling default.list")
  }

  if (!is.null(checkList)) {
    if (n1 != length(checkList)) stop("ERROR: in calling default.list")
    checkFlag <- 1
  } else {
    checkFlag <- 0
  } 

  if (is.null(inList)) inList <- list()

  listNames <- names(inList)
  for (i in 1:n1) {
    if (!(names[i] %in% listNames)) {
      if (!error[i]) {
        inList[[names[i]]] <- default[[i]]
      } else {
        temp <- paste("ERROR: the name ", names[i], " was not found", sep="")
        stop(temp)
      }
    } else if (checkFlag) {
      temp <- checkList[[i]]
      if (!all(is.na(temp))) {
        if (!all(inList[[names[i]]] %in% checkList[[i]])) {
          temp <- paste("ERROR: the name '", names[i], 
                      "' has an invalid value", sep="")
          stop(temp)
        }
      }
    }
  }

  inList

} # END: default.list

wrtTab <- function(x, f) {
  write.table(x, file=f, sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
}

addParmNames <- function(obj, retlist, retlistNames, op) {

  tmp <- retlistNames %in% names(retlist)
  retlistNames <- retlistNames[tmp]

  n <- length(retlistNames) 
  if (!n) return(retlist)

  # Get parm names
  parms     <- colnames(obj$z.all)
  test      <- obj$test.cols.2
  nuis      <- obj$nuis.cols.2
  ntest     <- length(test)
  nnuis     <- length(nuis)
  testParms <- parms[test]
  nuisParms <- parms[nuis]
  colNames <- c("score.estimates", "score.information")
  Names    <- c("stage2.estimates")
  fixed    <- op$fixed.effects
  
  for (name in retlistNames) {
    x <- retlist[[name, exact=TRUE]]
    if (name %in% colNames) {
      colnames(x) <- testParms
    } else if (name %in% Names) {
      if (fixed) {
        names(x) <- nuisParms
      } else {
        names(x) <- parms
      }
    }
    retlist[[name]] <- x
  }

  retlist

} # END: addParmNames
