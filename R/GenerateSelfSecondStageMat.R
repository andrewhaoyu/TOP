
#' Title
#'
#' @param z.design
#' @param x
#' @param M
#' @param full.second.stage.names
#' @param covar.names
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
GenerateSelfSecondStageMat <- function(x.self.design,
                                       z.design,
                                      M,
                                     full.second.stage.names,
                                     delta){

  x <- x.self.design
  ###1 for intercept
  if(is.vector(x)){
    x <- matrix(x,ncol=1)
    covar.names <- "Noname"
  }else{
    covar.names <- colnames(x)
  }
  total.covar.number <- 1+ncol(x)
  second.stage.cat <- ncol(z.design)

  max.number.second.stage.parameter <- second.stage.cat

  second.stage.mat <- matrix(NA,max.number.second.stage.parameter,(total.covar.number-1))

  delta.no.inter <- delta[(M+1):length(delta)]
  ind.delta <- 0
  ind.covar <- 0



    for(i in 1:(total.covar.number-1)){
      ind.covar <- ind.covar+1
      second.stage.mat[1:second.stage.cat,ind.covar]<-
        delta.no.inter[(ind.delta+1):(ind.delta+second.stage.cat)]
      ind.delta <- ind.delta + second.stage.cat
    }




  colnames(second.stage.mat) <- covar.names
  rownames(second.stage.mat) <- full.second.stage.names[1:max.number.second.stage.parameter]
  places = 3
  second.stage.mat <- round(second.stage.mat,digits = 3)
  return(second.stage.mat)

}
