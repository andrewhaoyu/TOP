
#' Title
#'
#' @param z.design.additive
#' @param M
#' @param tumor.names
#'
#' @return
#' @export
#'
#' @examples
GenerateSubtypesName <- function(z.design.additive,
                                 M,
                                 tumor.names) {
  z.design.standard <- z.design.additive[,2:ncol(z.design.additive)]
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



