

scan.score <- function(geno.list, pheno.list, scan.options=NULL) {

  geno.list    <- check.geno.list(geno.list)
  pheno.list   <- check.pheno.list(pheno.list)
  options      <- check.scan.op(scan.options, geno.list)
  scan.options <- NULL

  # Read in the phenotype data 
  x <- readPheno(pheno.list)

  # Get the genotype subject ids
  gids <- getIds.file(geno.list$subject.list) 
  geno.list$n.subjects <- length(gids)

  # Subset x based on genotype ids
  tmp <- x[, pheno.list$id.var] %in% gids
  x   <- x[tmp, , drop=FALSE]
  if (!nrow(x)) stop("ERROR: no matching ids in the phenotype and genotype data")

  # Check for errors in the data
  checkDataAndVars(x, pheno.list$cc.var, pheno.list$subtype.vars, NULL, pheno.list$covars.obj, check.test=FALSE)

  # Get all objects we need. All missing values will be removed
  obj0 <- getAllObjects(x, pheno.list$cc.var, pheno.list$subtype.vars, NULL, pheno.list$covars.obj, options)

  # Get the order matching the subject ids in phenotype and genotype data
  order0 <- match(rownames(obj0$y.obj), gids)
  tmp    <- rownames(x) %in% rownames(obj0$y.obj)
  x      <- x[tmp, , drop=FALSE]
  if (!nrow(x)) stop("ERROR: no subjects to include in the analysis")

  # Add a column for each SNP. This will  be needed if SNPs have missing values
  x    <- addSNP.col(x)

  options$snpv  <- colnames(x)[length(colnames(x))]
  delta0Flag    <- (!options$userDefined.delta0) && (!is.null(options[["delta0.list", exact=TRUE]]))
  if (delta0Flag) {
    options$delta0.list$test.cols.2 <- obj0$test.cols.2
    options$delta0.list$nuis.cols.2 <- obj0$nuis.cols.2
  }

  # Get the vector to output for each snp
  out.vec              <- getOutVec(options)
  n.out                <- length(out.vec)
  options$delta.nuis   <- obj0[["delta.nuis", exact=TRUE]]  
  obj0$delta.nuis      <- NULL

  index   <- geno.list$start.row - 1
  stopRow <- geno.list$stop.row

  # Open the genotype file and go to correct row
  fid <- getFilePtr(geno.list)

  # Open the output data
  ofid <- openOut(options$out.file, out.vec)
  
  while (1) {
    vec   <- scan(fid, what="character", sep=" ", nlines=1, quiet=TRUE)
    if (!length(vec)) break  # End of file
    index <- index + 1

    info  <- try(getSnpInfo(vec, geno.list, options))    
    if (checkTryError(info)) next  # Error with this row
    if (!info$use) next
    
    dosage  <- info$dosage[order0]
    tmp     <- try(scan_mixed(dosage, x, pheno.list, obj0, gids, options))
    eflag   <- checkTryError(tmp)
    out.vec <- setOutVec(out.vec, options, info, tmp)
    write(out.vec, file=ofid, ncolumns=n.out, sep="\t")

    if (index >= stopRow) break
    if (delta0Flag && !eflag && tmp$add.snp) {
      options$delta0.list <- update.delta0.list(dosage, tmp$stage2.estimates, options$delta0.list) 
    }

  } # END: while
  close(fid)

  out.vec[] <- "*"
  write(out.vec, file=ofid, ncolumns=n.out, sep="\t")

  if (!is.null(ofid)) close(ofid)

  NULL

} # END: scan.score

scan_mixed <- function(dosage, x, pheno.list, obj0, gids, options) {

  # dosage is already ordered
  if (is.null(ncol(dosage))) dim(dosage) <- c(length(dosage), 1)
  missFlag <- any(is.na(dosage))
  snpv     <- options$snpv
  add      <- FALSE
  flip     <- FALSE

  if (options$fixed) {
    if (missFlag) {
      x[, snpv] <- dosage
      obj1   <- getAllObjects(x, pheno.list$cc.var, pheno.list$subtype.vars, snpv,
                    pheno.list$covars.obj, options)
      ret    <- fixedEffect.main(dosage[obj1$subset, , drop=FALSE], obj1, options, add.parm.names=FALSE)
    } else {
      ret    <- fixedEffect.main(dosage, obj0, options, add.parm.names=FALSE)
    }
  } else {
    x[, snpv] <- dosage

    obj1   <- getAllObjects(x, pheno.list$cc.var, pheno.list$subtype.vars, snpv,
                  pheno.list$covars.obj, options)
    add    <- obj1$support.obj$add.snp
    flip   <- obj1$support.obj$flip.snp
    ret    <- randomEffect.main(dosage[obj1$subset, , drop=FALSE], obj1, options, add.parm.names=FALSE)
  }
  ret$add.snp  <- add
  ret$flip.snp <- flip

  ret

} # END: scan_mixed

setupDeltaList <- function(scan.op) {

  ret <- scan.op$delta0.list
  ret <- default.list(ret, 
                      c("save.N", "max.r2"), 
                      list(50, 0.5))
  ret$deltaN  <- 0
  ret$add.col <- 0
  if (ret$save.N < 0) stop("ERROR: save.N must be >= 0")
  if ((ret$max.r2 < 0) || (ret$max.r2 > 1)) stop("ERROR: max.r2 must be between 0 and 1")
  if (ret$save.N == 0) ret <- NULL
  
  ret

} # END: setupDeltaList

update.delta0.list <- function(newsnp, newdelta, deltaList) {

  deltaN <- deltaList$deltaN
  save.N <- deltaList$save.N

  if (!deltaN) {
    deltaMat      <- matrix(NA, nrow=length(newdelta), ncol=save.N)
    deltaMat[, 1] <- newdelta
    snpMat        <- matrix(NA, nrow=length(newsnp), ncol=save.N)
    snpMat[, 1]   <- newsnp
    add.col       <- 1
    deltaN        <- 1
  } else {
    snpMat   <- deltaList$snpMat
    deltaMat <- deltaList$deltaMat
    add.col  <- deltaList$add.col
    add.col  <- add.col + 1
    if (add.col > save.N) add.col <- 1  # Go back to beginning
    snpMat[, add.col]   <- newsnp
    deltaMat[, add.col] <- newdelta
    deltaN              <- min(deltaN + 1, save.N)
  }

  deltaList$deltaMat <- deltaMat
  deltaList$snpMat   <- snpMat
  deltaList$add.col  <- add.col
  deltaList$deltaN   <- deltaN

  deltaList

} # END: update.delta0.list


getDelta.scan <- function(newsnp, scan.op, obj0.delta0, nc) {

  # Check for a user defined delta0
  delta0 <- scan.op$options[["delta0", exact=TRUE]]
  flag   <- scan.op$userDefined.delta0
  if (is.null(nc)) nc <- 0
  del    <- NULL
  dlist  <- scan.op$delta0.list
  add    <- FALSE

  # We need to return a list, since getDelta.R2 does
  if (flag) {
    # Use the user-defined delta0
    del <- delta0
  } else if (dlist$deltaN > -1) {
    # Use LD info if it exists
    if (dlist$deltaN) {
      ret <- getDelta.R2(newsnp, dlist)
      add <- ret$add.snp.flag
      del <- ret$delta
    } else {
      add <- TRUE # First one must be added
    }
  } else {
    del <- obj0.delta0
  }
  len <- length(del)

  # Check that the length matches the z.all matrix
  if (nc && len && (nc != len)) {
    del2 <- rep(0, nc)
    if (len < nc) {
      # Recall that there are no fit parms for the test cols, so add in zeros
      nids <- dlist$nuis.cols.2
      if (len == length(nids)) del2[nids] <- del 
    } 
    del <- del2
  }

  list(delta=del, add.snp.flag=add)

} # END: getDelta.scan

addSNP.col <- function(x) {

  new      <- "_SNP_0g7h4d9h9hk47f"
  x[, new] <- NA
  x

} # END: addSNP.col

getSNP.info <- function(vec, geno.list) {

  ret <- list(use=FALSE, imputed=vec[1], snp=vec[2], loc=vec[3], 
              a1=vec[4], a2=vec[5], error="")
  ret

} # END: 

getExpDosage <- function(vec) {

  len   <- length(vec)
  ids   <- seq(from=1, to=len, by=3)
  svec1 <- vec[ids]
  svec2 <- vec[ids+1]
  svec3 <- vec[ids+2]

  # Check for missing values
  temp <- svec1 + svec2 + svec3 == 0
  temp[is.na(temp)] <- TRUE
  if (any(temp)) {
    svec1[temp] <- NA
    svec2[temp] <- NA
    svec3[temp] <- NA
  }

  svec2 + 2*svec3

} # END: getExpDosage

checkNextSNP <- function(snp.info, geno.list, scan.options) {

  min.maf <- scan.options$min.maf
  if (min.maf) {
    eaf <- snp.info$eaf
    maf <- min(eaf, 1-eaf)
    if (maf < min.maf) snp.info$use=FALSE
  }

  snp.info

} # END: checkNextSNP

getSnpInfo <- function(vec, geno.list, scan.options) {

  n   <- length(vec)
  if (!n) return(NULL)
  if (n < 5) stop("ERROR with genotype row")
  ret   <- getSNP.info(vec, geno.list)
  nsub  <- geno.list$n.subjects
  if ((n-5)/3 != nsub) stop("ERROR: Incorrect number of genotypes")
  
  vec        <- as.numeric(vec[-(1:5)])
  dosage     <- getExpDosage(vec)
  ret$eaf    <- 0.5*mean(dosage, na.rm=TRUE)
  ret$dosage <- dosage
  ret$use    <- TRUE
  ret        <- checkNextSNP(ret, geno.list, scan.options)

  ret

} # END: getNextSNP

openOut <- function(f, out.vec) {

  ret <- file(f, "w")
  write(names(out.vec), file=ret, ncolumns=length(out.vec), sep="\t")
  ret

} # END: openOut

getOutVec <- function(scan.op) {

  cols <- c("SNP", "Position", "RefAllele", "EffectAllele", "EAF")
  cols <- c(cols, "Pvalue.GTA.fixed")
  if (!scan.op$fixed) cols <- c(cols, "Pvalue.GTA.mixed")
  cols <- c(cols, "Error")
  ret <- rep("", length(cols))
  names(ret) <- cols
  ret

} # END: getOutVec

setOutVec <- function(out.vec, scan.op, info, res) {

  ecol          <- length(out.vec)
  out.vec[]     <- ""
  out.vec[1:5]  <- c(info$snp, info$loc, info$a1, info$a2, info$eaf)
  out.vec[ecol] <- info$error

  if (is.list(res)) {
    scr <- res[["global.score.test", exact=TRUE]]
    if (!is.null(scr)) {
      if (scan.op$fixed) {
        out.vec["Pvalue.GTA.fixed"] <- scr$p.value
      } else {
        out.vec["Pvalue.GTA.mixed"] <- scr$p.value
        out.vec["Pvalue.GTA.fixed"] <- scr$p.value.fixed
        msg <- scr[["error.message", exact=TRUE]]
        if (!is.null(msg)) out.vec[ecol] <- msg
      }
    }
  } else {
    if (is.null(res)) {
      out.vec[ecol] <- "Error occurred"
    } else if ((is.character(res)) && (length(res) == 1)) {
      res <- gsub("\n", " ", res, fixed=TRUE)
      res <- gsub("\t", " ", res, fixed=TRUE)
      res <- gsub("\r", " ", res, fixed=TRUE)
      out.vec[ecol] <- res
    }
  }

  out.vec

} # END: setOutVec

readPheno <- function(pheno.list) {
 
  x <- read.table(pheno.list$file, header=pheno.list$header, 
                  sep=pheno.list$delimiter, as.is=TRUE)

  tmp <- pheno.list[["subsetData", exact=TRUE]]
  if (!is.null(tmp)) x <- subsetData.list(x, tmp)

  # Keep only the variables we need
  vars <- pheno.list[["keep.vars", exact=TRUE]]
  if (!is.null(vars)) x <- x[, vars, drop=FALSE]
  x    <- as.data.frame(x, stringsAsFactors=FALSE)
  
  # Assign the subject ids as row names
  pids <- x[, pheno.list$id.var]
  if (any(duplicated(pids))) stop("ERROR: pheno.list$id.var must have all unique ids")
  rownames(x) <- pids

  x

} # END: readPheno

check.subject.list <- function(slist, format=NULL) {

  len   <- length(slist)
  if (!len) return(slist)

  if ((len == 1) && (is.character(slist))) {
    # Assume a file
    slist <- list(file=slist) 
  }
  if (!is.list(slist)) return(slist)
  
  ff <- slist[["file", exact=TRUE]]  
  if (is.null(ff)) return(slist)

  # Check that the file exists
  if (!file.exists(ff)) stop("ERROR with subject.list")

  header <- slist[["header", exact=TRUE]]
  if (is.null(header)) slist$header <- 0

  sep  <- slist[["delimiter", exact=TRUE]]
  if (is.null(sep)) slist$delimiter <- ""

  # If id.var is specified return
  id.var <- slist[["id.var", exact=TRUE]]
  if (!is.null(id.var)) return(slist)

  # Get the number of columns
  row1 <- scan(ff, what="character", sep="", nlines=1, quiet=TRUE)
  nc   <- length(row1)

  if (nc == 1) {
    slist$id.var <- -1
  } else {
    slist$id.var <- 1
  }
  
  slist

} # END: check.subject.list

check.geno.list <- function(geno.list) {

  if (is.null(geno.list)) stop("ERROR: geno.list must be specified")
  if (!is.list(geno.list)) stop("ERROR: geno.list must be a list")

  f <- geno.list[["file", exact=TRUE]]
  if (is.null(f)) stop("ERROR: geno.list$file must be specified")
  sublist <- geno.list[["subject.list", exact=TRUE]]
  if (is.null(sublist)) stop("ERROR: geno.list$subject.list must be specified")
  f <- sublist[["file", exact=TRUE]]
  if (is.null(f)) stop("ERROR: subject.list$file must be specified")
  geno.list$subject.list <- check.subject.list(sublist)
  
  geno.list <- default.list(geno.list, c("start.row", "stop.row"),
                     list(1, Inf), error=c(0, 0))
  if (geno.list$stop.row < 1) geno.list$stop.row <- Inf
  if (geno.list$start.row < 1) stop("ERROR: geno.list$start.row must be >= 1")
  if (geno.list$stop.row < geno.list$start.row) stop("ERROR: geno.list$stop.row must be >= geno.list$start.row")
  
  geno.list$format <- "impute"

  geno.list

} # END: check.geno.list

check.vars.file <- function(f, sep, vars) {

  miss  <- NULL
  vec   <- scan(f, what="character", nlines=1, sep=sep, quiet=TRUE)
  tmp   <- !(vars %in% vec)
  if (any(tmp)) miss <- vars[tmp]
  
  miss

} # END: chck.vars.file

check.pheno.list <- function(pheno.list) {

  if (is.null(pheno.list)) stop("ERROR: pheno.list must be specified")
  if (!is.list(pheno.list)) stop("ERROR: pheno.list must be a list")
  f <- pheno.list[["file", exact=TRUE]]
  if (is.null(f))  stop("ERROR: pheno.list$file must be specified")
  if (!file.exists(f)) {
    str <- paste("ERROR: the file ", f, " does not exist.", sep="")
    stop(str)
  }
  if (is.null(pheno.list[["id.var", exact=TRUE]])) stop("ERROR: pheno.list$id.var must be specified")
  if (is.null(pheno.list[["cc.var", exact=TRUE]])) stop("ERROR: pheno.list$cc.var must be specified")
  if (is.null(pheno.list[["subtype.vars", exact=TRUE]])) stop("ERROR: pheno.list$subtype.vars must be specified")
  pheno.list$subtype.vars <- unique(pheno.list$subtype.vars)
  pheno.list <- default.list(pheno.list, c("delimiter", "header"), list("", 1))

  tmp        <- getAllVars(pheno.list)
  vars       <- tmp$all.vars
  pheno.list$all.vars  <- vars
  pheno.list$keep.vars <- tmp$keep.vars

  miss <- check.vars.file(f, pheno.list$delimiter, vars)
  if (length(miss)) {
    str <- paste(miss, collapse=", ", sep="")
    msg <- paste("ERROR: the variables ", str, " were not found in pheno.list$file", sep="")
    stop(msg)
  }

  slist <- pheno.list[["subsetData", exact=TRUE]]
  if (!is.null(slist)) {
    slist <- check.subsetData.list(slist) 
  } 
  
  pheno.list

} # END: check.pheno.list

check.scan.op <- function(scan.op, geno.list) {

  # Return all options in one list
  scan.op <- default.list(scan.op, 
                          c("method", "min.maf"), 
                          list("fixed", 0), 
                          error=c(0, 0))
  if (!(scan.op$method %in% c("fixed", "random"))) stop("ERROR: scan.options$method must be fixed or random")

  fixed                      <- scan.op$method == "fixed"
  options                    <- scan.op[["options", exact=TRUE]]
  delta0                     <- options[["delta0", exact=TRUE]]
  options                    <- check_options(options, fixed=fixed)
  options$stopOnDelta0.error <- FALSE
  options$method             <- scan.op$method
  options$min.maf            <- scan.op$min.maf

  out <- scan.op[["out.file", exact=TRUE]]
  if (is.null(out)) {
    a   <- geno.list$start.row
    b   <- geno.list$stop.row
    out <- paste(getwd(), "/", basename(geno.list$file), "_out_", sep="")
    if ((a > 1) || !is.finite(b)) {
      if (!is.finite(b)) b <- "Inf"
      out <- paste(out, a, "_", b, sep="")
    }
  }
  options$out.file           <- out
  options$delta0.list        <- setupDeltaList(scan.op)
  options$userDefined.delta0 <- !is.null(delta0)

  options

} # END: check.scan.op

# Function to return the variables from a particular object
getAllVars <- function(obj, names=c("cc.var", "subtype.vars", "covars.obj", "id.var")) {

  if (is.null(obj)) return(NULL)
  clss <- class(obj)
  keep <- NULL

  # Character vector
  if ((is.vector(obj)) && ("character" %in% clss)) {
    vars <- unique(obj)
    ret  <- list(all.vars=vars, keep.vars=vars)
    return(ret)
  }

  # Formula
  if ("formula" %in% clss) {
    vars <- unique(all.vars(obj))
    ret  <- list(all.vars=vars, keep.vars=vars)
    return(ret)
  }

  # List
  if ("list" %in% clss) {
    ret <- NULL
    if (is.null(names)) names <- 1:length(obj)
    for (nn in names) {
      obj.n <- obj[[nn, exact=TRUE]]
      if (is.null(obj.n)) next

      clss <- class(obj.n)      
      if ((is.vector(obj)) && ("character" %in% clss)) {
        ret <- c(ret, obj.n)
      } else if ("formula" %in% clss) {
        ret <- c(ret, all.vars(obj.n))
      }
    }
    keep <- unique(ret)

    slist <- obj[["subsetData", exact=TRUE]]
    if (!is.null(slist)) {
      svars <- getSubsetDataVars(slist)
      ret   <- c(ret, svars)
    }
   
    return(list(all.vars=unique(ret), keep.vars=keep))
  }

  stop("ERROR in formulaVars: obj is of wrong type")

} # END: getAllVars

# Function to return the variables in a list of type subsetData
getSubsetDataVars <- function(slist) {

  n   <- length(slist)
  ret <- character(n)
  for (i in 1:n) {
    ret[i] <- slist[[i]]$var
  }

  ret

} # END: getSubsetDataVars

# Function to check a list of type subsetData
check.subsetData.list <- function(slist) {

  n <- length(slist)
  for (i in 1:n) {
    temp <- default.list(slist[[i]], c("var", "operator", "value"), 
                         list("ERROR", "ERROR", "ERROR"), error=c(1, 1, 1))
  }

  slist

} # END: check.subsetData.list

# Function to open file and postion file pointer
getFilePtr <- function(geno.list) {

  row <- geno.list$start.row
  fid <- gzfile(geno.list$file, "r")
  if (row > 1) {
    temp <- scan(fid, what="character", nlines=1, sep="\n", skip=row-2, quiet=TRUE)  
  }

  fid

} # END: getFilePtr

# Function to get a column of ids from a file
getIds.file <- function(file.list, rm.space=1, lower=0, upper=0, id.sep=":") {

  len  <- length(file.list)
  flag <- FALSE
  if (len == 1) {
    flag <- try(file.exists(file.list), silent=TRUE)
    if (checkTryError(flag)) flag <- FALSE
  } 

  if ((!flag) && (is.vector(file.list)) && (!("list" %in% class(file.list)))) {
    ids <- file.list
  } else {
    #file.list <- check.file.list(file.list)
    f         <- file.list$file
    sep       <- file.list$delimiter

    # Get the number of columns
    row1 <- scan(f, what="character", sep=sep, quiet=TRUE, nlines=1)
    nc   <- length(row1)
    
    if ((flag) && (nc == 1)) {
      file.list$id.var <- -1
      file.list$header <- 0
    } 
    if (is.null(file.list[["id.var", exact=TRUE]])) file.list$id.var <- 1
    id.var <- file.list$id.var
    if ((length(id.var) == 1) && (id.var %in% -1)) {
      ids <- scan(f, what="character", sep=sep, quiet=TRUE)
    } else {  
      # Allow for 2 id variables
      x   <- matrix(scan(f, what="character", sep=sep, quiet=TRUE), byrow=TRUE, ncol=nc)
      if (is.numeric(id.var)) {
        ids <- x[, id.var]
      } else {
        temp <- row1 %in% id.var
        if (sum(temp) != 1) stop("ERROR in getIds.file: with id.var")
        col <- (1:nc)[temp] 
        ids <- x[, col]
      }
      m <- ncol(ids)
      if (is.null(m)) m <- 0
      if (m == 2) ids <- paste(makeVector(ids[, 1]), id.sep, makeVector(ids[, 2]), sep="") 

      ids <- makeVector(ids) 
      if (file.list$header) ids <- ids[-1]
    }
  }
  if (rm.space) ids <- removeWhiteSpace(ids)
  if (lower) ids <- tolower(ids)
  if (upper) ids <- toupper(ids)

  ids

} # END: getIds.file

checkTryError <- function(obj) {

  ret <- 0
  classObj <- class(obj)
  if ("try-error" %in% classObj) ret <- 1

  ret

} # END: checkTryError

makeVector <- function(x) {

  d <- dim(x)
  if (is.null(d)) return(x)

  nn <- NULL
  if (d[1] == 1) {
    nn <- colnames(x)
  } else if (d[2] == 1) {
    nn <- rownames(x)
  }
  dim(x) <- NULL
  if (!is.null(nn)) names(x) <- nn
  if ((!is.vector(x)) && (is.list(x))) x <- unlist(x)
  x
 
} # END: makeVector

removeWhiteSpace <- function(str, leading=1, trailing=1) {

  if ((leading) && (trailing)) {
    ret <- gsub("^\\s+|\\s+$", "", str, perl=TRUE)
  } else if (leading) {
    ret <- gsub("^\\s+", "", str, perl=TRUE)
  } else if (trailing) {
    ret <- gsub("\\s+$", "", str, perl=TRUE)
  } else {
    ret <- str
  }

  ret

} # END: removeWhiteSpace

# Function to subset data by a single variable
subsetData.var <- function(data, var, operator, value, which=1, 
                           returnRows=0, na.value=FALSE) {

  # data
  # var
  # operator
  # value
  # which        1 or -1 to keep or drop
  #              The default is 1
  # returnRows   0 or 1
  #              The default is 0
  # na.value     TRUE or FALSE on what to do with NAs
  #              The default is FALSE

  lenv <- length(value)
  if (any(is.na(value))) {
    if (lenv > 1) stop("ERROR: IN subsetData.var: with NA in value")
    naFlag  <- 1
    numFlag <- 0
  } else {
    numFlag <- is.numeric(value)
    naFlag  <- 0
  }

  # Be careful with a vector for value
  if (lenv > 1) {
    if (numFlag) {
      value <- paste(value, collapse=",", sep="")
      value <- paste("c(", value, ")", sep="")
    } else {
      # Change the operator
      if (operator == "%in%") {
        operator <- "=="
      } else { 
        stop("ERROR: value cannot be a character vector")
      }
    }
  }

  if (naFlag) {
    temp <- ""
    if (operator == "!=") temp <- "!"  
    callStr <- paste("(", temp, "is.na(data[, var]))", sep=" ") 
  }
  else if (numFlag) {
    callStr <- paste("(as.numeric(data[, var])", operator, value, ")", sep=" ") 
  } else {
    callStr <- ""
    for (i in 1:lenv) {
      temp <- paste('(data[, var] ', operator, ' "', value[i], '")', sep='')
      if (i < lenv) temp <- paste(temp, " | ", sep="")
      callStr <- paste(callStr, temp, sep="")
    }
  }

  rows <- eval(parse(text=callStr))
  rows[is.na(rows)] <- na.value
  if (returnRows) return(rows)

  data <- removeOrKeepRows(data, rows, which=which)
  data

} # END: subsetData.var

# Function to subset data. Returns the subsetted data, or if op$which = 0,
#   returns the vector of TRUE/FALSE to subset
subsetData.list <- function(data, slist, returnRows=0) {

  # data
  ###################################################################
  # slist       List of sublists with names:
  #  var
  #  operator
  #  value
  #  which      1 or -1 to include or NOT within each list
  #             The default is 1
  #  logic.op   Logical operator 
  #             Only used starting from the second sublist
  #             The default is "&"
  #  na.value   TRUE or FALSE
  #             The default is FALSE
  #  last.which 1 or -1 If -1, the rows defined by all the sublists will
  #             be negated. The value for last.which in the last sublist
  #             is the only one applied.
  #             The default is 1.
  ####################################################################
  # returnRows  Set to 1 to return the logical vector of rows instead of
  #             the data.
  #             The default is 0

  n      <- length(slist)
  wvec   <- rep(9999, times=n)
  cnames <- colnames(data)
  cflag  <- !is.null(cnames)

  for (i in 1:n) {

    tlist    <- slist[[i]]
    if (!is.list(tlist)) stop("ERROR in subsetData.list: Input list is incorrect")
    tlist    <- default.list(tlist, 
                  c("var", "operator", "value", "which", "logic.op", "na.value", "last.which"), 
                  list("ERROR", "ERROR", "ERROR", 1, "&", FALSE, 1), error=c(1,1,1,0,0,0,0))
    var      <- tlist[["var", exact=TRUE]]
    if ((cflag) && (is.character(var))) {
      if (!(var %in% cnames)) {
        temp <- paste("ERROR in subsetData.list: ", var, " not in data", sep="")
        print(temp)
        stop()
      }
    }
    operator <- tlist[["operator", exact=TRUE]]
    value    <- tlist[["value", exact=TRUE]]
    which    <- tlist[["which", exact=TRUE]]
    last     <- tlist[["last.which", exact=TRUE]]
    wvec[i]  <- which
    na.value <- tlist[["na.value", exact=TRUE]]
    temp     <- subsetData.var(data, var, operator, value, returnRows=1, na.value=na.value)
    if (which == -1) temp <- !temp
    if (i == 1) {
      rows <- temp
    } else {
      logic.op <- tlist[["logic.op", exact=TRUE]]
      callStr  <- paste("rows ", logic.op, " temp", sep="") 
      rows     <- eval(parse(text=callStr))      
    }
  }

  if (last == -1) rows <- !rows
  if (returnRows) return(rows)

  data <- removeOrKeepRows(data, rows, which=1)

  data

} # END: subsetData.list

# Function to remove columns or variables from a matrix or data frame
removeOrKeepRows <- function(x, rows, which=1) {

  # x      Matrix or data frame 
  # rows   Vector of row numbers, character vector of row names,
  #        or logical vector.
  # which  1 or -1 to keep or remove rows
  #        The default is 1 

  n      <- dim(x)
  dfFlag <- is.data.frame(x)
  if (is.null(n)) stop("ERROR: x should be 2 dimensional")
  if (which != 1) which <- -1
  rnames <- rownames(x)

  if (is.logical(rows)) {
    if (length(rows) != n[1]) stop("ERROR with logical vector rows")
    if (which != 1) rows <- !rows
  } else if (is.character(rows)) {
    rows <- match(rows, rnames)
    if (any(is.na(rows))) stop("ERROR in removeOrKeepRows")
    rows <- which*rows
  } else {
    rows <- which*rows
  }
  cnames <- colnames(x)

  # Keep or remove rows
  x <- x[rows, ]
  
  # Check for NULL dimension
  if (is.null(dim(x))) {
    dim(x) <- c(length(x)/n[2], n[2])
    if (dfFlag) x <- data.frame(x)
    if (!is.null(cnames)) colnames(x) <- cnames
    rownames(x) <- rnames[rows]
  }

  x

} # END: removeOrKeepRows
