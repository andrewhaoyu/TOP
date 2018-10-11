#' Title
#'
#' @param A
#' @param B
#'
#' @return
#' @export
#'
#' @examples
rowmatch <- function(A,B) {
  # Rows in A that match the rows in B
  f <- function(...) paste(..., sep=":")
  if(!is.matrix(B)) B <- matrix(B, 1, length(B))
  a <- do.call("f", as.data.frame(A))
  b <- do.call("f", as.data.frame(B))
 which(a%in%b)
}



CleanNotEnoughSubtype <- function(y,
baselineonly=NULL,
additive=NULL,
pairwise.interaction=NULL,
saturated=NULL,
missingTumorIndicator = 888,
delta0= NULL){
  missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator)
  y.pheno.complete <- y[-missing.data.vec,]
  freq.subtypes <- GenerateFreqTable(y.pheno.complete)
  cutoff <- 10
  freq.table.mis <- GenerateFreqTable(y)
  freq.table.com <- GenerateFreqTable(y.pheno.complete)
  idx <- which(freq.table.com[,ncol(freq.table.com)]>cutoff)
  freq.table.com <- freq.table.com[idx,-ncol(freq.table.com)]
  n.col <- ncol(freq.table.mis)


  subtype.result <- NULL
  for(i in 1:nrow(freq.table.mis)){
    record <- 0
    temp <- freq.table.mis[i,1:(n.col-1),drop=F]
    idx <- which(temp==missingTumorIndicator)

    if(length(idx)!=0){
      temp <- temp[-idx,drop=F]
      freq.table.ref <- freq.table.com[,-idx,drop=F]
      for(j in 1:nrow(freq.table.ref)){
        if(all.equal(as.vector(temp),as.vector(freq.table.ref[j,]))==1){
          record <- 1
        }
      }
      if(record == 0){
        subtype.result <- rbind(subtype.result,freq.table.mis[i,1:(n.col-1)])
      }

    }else{

      freq.table.ref <- freq.table.com
      for(j in 1:nrow(freq.table.ref)){
        if(all.equal(as.vector(temp),as.vector(freq.table.ref[j,]))==1){
          record <- 1
        }
      }
      if(record == 0){
        subtype.result <- rbind(subtype.result,freq.table.mis[i,1:(n.col-1)])
      }

    }


  }



}
