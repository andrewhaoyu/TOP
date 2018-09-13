

#' Title
#'
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param x.self.design
#'
#' @return
#' @export
#'
#' @examples
GenerateCovarName <- function(baselineonly,
                              additive,
                              pairwise.interaction,
                              saturated,
                              x.self.design=NULL){

  full.names <- NULL

  if(is.null(x.self.design)==0){
    if(is.vector(x.self.design)==1){
      x.self.design = as.matrix(x.self.design)
      colnames(x.self.design) = "self design variable"
      full.names<- c(full.names,colnames(x.self.design))
    }else{
      full.names<- c(full.names,colnames(x.self.design))
    }

  }

  if(is.null(baselineonly)==0){
    if(is.vector(baselineonly)==1){
      baselineonly = as.matrix(baselineonly)
      colnames(baselineonly) = "baselineonlly variable"
      full.names<- c(full.names,colnames(baselineonly))
    }else{
      full.names<- c(full.names,colnames(baselineonly))
    }

  }
  if(is.null(additive)==0){

    if(is.vector(additive)==1){
      additive = as.matrix(additive)
      colnames(additive) = "additive variable"
      full.names<- c(full.names,colnames(additive))
    }else{
      full.names<- c(full.names,colnames(additive))
    }
  }
  if(is.null(pairwise.interaction)==0){

    if(is.vector(pairwise.interaction)==1){
      pairwise.interaction = as.matrix(pairwise.interaction)
      colnames(pairwise.interaction) = "pairwise.interaction variable"
      full.names<- c(full.names,colnames(pairwise.interaction))
    }else{
      full.names<- c(full.names,colnames(pairwise.interaction))
    }
  }
  if(is.null(saturated)==0){
    if(is.vector(saturated)==1){
      saturated = as.matrix(saturated)
      colnames(saturated) = "saturated variable"
      full.names<- c(full.names,colnames(saturated))
    }else{
      full.names<- c(full.names,colnames(saturated))
    }
  }
  return(full.names)
}
