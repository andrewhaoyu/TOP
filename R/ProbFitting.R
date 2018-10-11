####Calculate the conditional probability for E step
####output: missing.vec indicates the subject index with missing tumor
####output: missing.mat is a binary matrix indicating the potential subtypes this missing person could be
####output: y_em the initial results for E step matrix
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
ProbFitting <- function(delta0,y,x.all,z.standard,z.all,missingTumorIndicator=NULL){
  if(is.null(missingTumorIndicator)==1){
    n <- nrow(y)
    M <- nrow(z.standard)
    y_em <- matrix(0,nrow=n,ncol = M)
    beta <- matrix(z.all%*%delta0,ncol = M)
    ##add intercept to x.all
    x.all.inter <- as.matrix(cbind(1,x.all))
    index = 1
    for(i in 1:nrow(y)){
      if(y[i,1]==1) {
        ###find out which tumor characteristic is observed
        ###-3.14 is just a random number to make the algorithm run
        ###since there is no missing, all of the 2:ncol(y) will be chose
        idx <- 1:(ncol(y)-1)
        ###jdx is the potential subtype this missing person could be
        jdx <- apply(z.standard,1,function(t){all(t[idx]==y[i,idx+1])})
        jdx <- which(jdx==T)
        ####get the conditional probability
        y_em[i,jdx] <- 1
      }
    }
    return(list(y_em=y_em))
  }else{
    n <- nrow(y)
    M <- nrow(z.standard)
    y_em <- matrix(0,nrow=n,ncol = M)
    missing.vec = rep(0,n)
    missing.mat = matrix(0,nrow=n,ncol=M)

    beta <- matrix(z.all%*%delta0,ncol = M)
    ##add intercept to x.all
    x.all.inter <- as.matrix(cbind(1,x.all))
    index = 1


    for(i in 1:nrow(y)){
      if(y[i,1]==1) {
        ###find out which tumor characteristic is observed
        idx <- which(y[i,2:ncol(y)]!=missingTumorIndicator)
        ###jdx is the potential subtype this missing person could be
        jdx <- apply(z.standard,1,function(t){all(t[idx]==y[i,idx+1])})
        if(sum(jdx)!=1){
          missing.vec[index] <- i
          missing.mat[index,] <- jdx
          index <- index + 1
          jdx <- which(jdx==T)
          ####get the conditional probability
          temp <- exp(x.all.inter[i,]%*%beta[,jdx])
          y_em[i,jdx] <- temp/sum(temp)
        }else{
          y_em[i,jdx] <- 1
        }

      }
    }

    missing.vec <- missing.vec[1:(index-1)]
    missing.mat <- missing.mat[1:(index-1),]
    complete.vec <- c(1:n)
    complete.vec <- complete.vec[!(complete.vec%in%missing.vec)]

    return(list(y_em=y_em,missing.vec = missing.vec , missing.mat = missing.mat,complete.vec = complete.vec))
  }

}





# Loglikelihood <- function(delta,y,x.all,z.standard,z.all,missingTumorIndicator){
#   n <- nrow(y)
#   M <- nrow(z.standard)
#   p.o <- rep(0,n)
#   beta <- matrix(z.all%*%delta,ncol = M)
#   ##add intercept to x.all
#   x.all.inter <- as.matrix(cbind(1,x.all))
#   index = 1
#   for(i in 1:nrow(y)){
#     ###p.all represent 1+sum(exp(x*beta))
#     p.all <- 1
#     p.vec <- rep(0,M)
#     for(j in 1:M){
#       p.vec [j]<- exp(x.all.inter[i,]%*%beta[,j])
#       p.all <- p.all+p.vec[j]
#     }
#
#     if(y[i,1]==1) {
#       ###find out which tumor characteristic is observed
#       idx <- which(y[i,2:ncol(y)]!=missingTumorIndicator)
#       ###jdx is the potential subtype this missing person could be
#       jdx <- apply(z.standard,1,function(t){all(t[idx]==y[i,idx+1])})
#       temp <- t(p.vec)%*%jdx
#      p.o[i] <- temp/(p.all)
#     }else{
#      p.o[i] <- 1/p.all
#     }
#   }
#
#
#   return(p.o)
# }
