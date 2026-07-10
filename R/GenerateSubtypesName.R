
#' Title
#'
#' @param z.design.additive Input for `z.design.additive`.
#' @param M Input for `M`.
#' @param tumor.names Input for `tumor.names`.
#'
#' @return A result from this computation.
#' @export
#'
GenerateSubtypesName <- function(z.design.additive,
                                 M,
                                 tumor.names) {
  z.design.standard <- z.design.additive[,2:ncol(z.design.additive),drop=F]
  M <- nrow(z.design.standard)
  subtypes.names <- rep("",M)
  for(i in 1:M){
    for(j in 1:length(tumor.names)){
      temp <- paste0(tumor.names[j],z.design.standard[i,j])
      subtypes.names[i] <- paste0(subtypes.names[i],temp)
    }

  }
  return(subtypes.names)
}



