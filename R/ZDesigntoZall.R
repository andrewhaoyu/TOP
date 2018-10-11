#' Title
#'
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param z.design.baselineonly
#' @param z.design.additive
#' @param z.design.pairwise.interaction
#' @param z.design.saturated
#'
#' @return
#' @export
#'
#' @examples
ZDesigntoZall <- function(baselineonly,
                          additive,
                          pairwise.interaction,
                          saturated,
                          z.design.baselineonly,
                          z.design.additive,
                          z.design.pairwise.interaction,
                          z.design.saturated) {
  M <- nrow(z.design.additive)
  ##the number of covariates in different potential model structures
  baselineonly.number <- CountCovarNumber(baselineonly)
  additive.number <- CountCovarNumber(additive)
  pairwise.interaction.number <- CountCovarNumber(pairwise.interaction)
  saturated.number <- CountCovarNumber(saturated)
  ###second.stage.category for different model structures
  baselineonly.second.cat <- 1
  additive.second.cat <- ncol(z.design.additive)
  pairwise.interaction.second.cat <- ncol(z.design.pairwise.interaction)
  saturated.second.cat <- ncol(z.design.saturated)
  ###1 for intercept
  total.covar.number <- 1+ baselineonly.number+additive.number+
    pairwise.interaction.number+saturated.number

  z.all <- matrix(0,nrow=(M*total.covar.number),ncol = (M+baselineonly.number*baselineonly.second.cat+
                                                          additive.second.cat*additive.number+
                                                          pairwise.interaction.second.cat*pairwise.interaction.number)+saturated.second.cat*saturated.number)

  for(i in c("intercept","baselineonly",
             "additive","pairwise.interaction",
             "satuared")){
    ##we always keep intercept as saturated model and to simply, we always use diagnonal matrix for intercept
    if(i=="intercept"){
      ###row start and column start point for this category
      row.start <- 0
      column.start <- 0
      for(j in 1:M){
        z.all[row.start+1+(j-1)*total.covar.number,(column.start+j)] = 1
      }
    }else if(i=="baselineonly"){
      column.start = M
      row.start <- 1
      ###test whether there is any baselineonly variable
      if(baselineonly.number!=0){
        for(j in 1:M){
          for(k in 1:baselineonly.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*baselineonly.second.cat+1):
                    (column.start+k*baselineonly.second.cat)] <- as.matrix(z.design.baselineonly[j,])
          }
        }
      }
    }else if(i=="additive"){
      column.start <- M+baselineonly.number
      row.start <- 1+baselineonly.number
      if(additive.number!=0){
        for(j in 1:M){
          for(k in 1:additive.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*additive.second.cat+1):
                    (column.start+k*additive.second.cat)] <- as.matrix(z.design.additive[j,])
          }
        }
      }
    }else if(i == "pairwise.interaction"){
      column.start <- M+baselineonly.number+additive.number*additive.second.cat
      row.start <- 1+baselineonly.number+additive.number
      if(pairwise.interaction.number!=0){
        for(j in 1:M){
          for(k in 1:pairwise.interaction.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*pairwise.interaction.second.cat+1):
                    (column.start+k*pairwise.interaction.second.cat)] <- as.matrix(z.design.pairwise.interaction[j,])
          }
        }
      }
    }else {
      column.start <- M+baselineonly.number+additive.number*additive.second.cat+
        pairwise.interaction.number*pairwise.interaction.second.cat
      row.start <- 1+baselineonly.number+additive.number+pairwise.interaction.number
      if(saturated.number!=0){
        for(j in 1:M){
          for(k in 1:saturated.number){
            z.all[row.start+k+(j-1)*total.covar.number,
                  (column.start+(k-1)*saturated.second.cat+1):
                    (column.start+k*saturated.second.cat)] <- as.matrix(z.design.saturated[j,])
          }
        }
      }
    }



  }
  return(z.all)
}
