###Generate the potential Z Design Matrix for the four potential models


#' Title
#'
#' @param tumor.character.cat
#' @param tumor.number
#' @param tumor.names
#' @param freq.subtypes
#'
#' @return
#' @export
#'
#' @examples
GenerateZDesignBaselineonly <- function(tumor.character.cat,tumor.number,tumor.names,freq.subtypes){
  M = 1
  cutoff <- 10
  for(i in 1:tumor.number){
    M = M*length(tumor.character.cat[[i]])
  }

  z.design.baselineonly <- matrix(rep(1,M),M,1)

  freq = freq.subtypes[,ncol(freq.subtypes)]
  idx <- which(freq<=cutoff)
  if(length(idx!=0)){
    return(z.design.baselineonly[-idx,,drop=F])
  }else{
    return(z.design.baselineonly)
  }

}

###Generate z design matrix for main effect model

#' Title
#'
#' @param tumor.character.cat
#' @param tumor.number
#' @param tumor.names
#' @param freq.subtypes
#'
#' @return
#' @export
#'
#' @examples
GenerateZDesignAdditive <- function(tumor.character.cat,tumor.number,tumor.names,freq.subtypes){
  z.design.additive.text <- NULL
  cutoff <- 10
  for(i in 1:tumor.number){
    if(i==tumor.number){
      z.design.additive.text <- paste0(z.design.additive.text,
                                          "tumor.character.cat[[",i,"]]")
    }else{
      z.design.additive.text <- paste0(z.design.additive.text,
                                          "tumor.character.cat[[",i,"]],")
    }
  }
  z.design.additive.text <- paste0("z.design.additive <- expand.grid(",
                                      z.design.additive.text,
                                      ")")
  eval(parse(text=z.design.additive.text))
  z.design.additive <- cbind(1,z.design.additive)
  colnames(z.design.additive) <- GenerateZDesignNamesAdditive(tumor.names)
  freq = freq.subtypes[,ncol(freq.subtypes)]
  idx <- which(freq<=cutoff)
  if(length(idx!=0)){
    return(z.design.additive[-idx,])
  }else{
    return(z.design.additive)
  }
}

###Generate tumor.names for z design matrix for main effect

#' Title
#'
#' @param tumor.names
#'
#' @return
#' @export
#'
#' @examples
GenerateZDesignNamesAdditive <- function(tumor.names){

  z.design.names.additive <- "baseline effect"

  z.design.names.additive <- c(z.design.names.additive,
                                  paste0(tumor.names," main effect"))

  return(z.design.names.additive)

}

###Generate z design matrix for pairwise interaction model

#' Title
#'
#' @param tumor.character.cat
#' @param tumor.number
#' @param tumor.names
#' @param freq.subtypes
#'
#' @return
#' @export
#'
#' @examples
GenerateZDesignPairwiseInteraction <- function(tumor.character.cat,tumor.number,tumor.names,freq.subtypes){
  cutoff <- 10
  z.design.pairwise.interaction <-
    GenerateZDesignAdditive(tumor.character.cat,
                              tumor.number,
                              tumor.names,
                              freq.subtypes)
  z.design.names.pairwise.interaction <- colnames(z.design.pairwise.interaction)
  all.pairwise.combnation <- combn(tumor.number,2)+1
  combn.number <- ncol(all.pairwise.combnation)

  for(i in 1:combn.number){
    col1 <- all.pairwise.combnation[1,i]
    col2 <- all.pairwise.combnation[2,i]
    newcol <-  z.design.pairwise.interaction[,col1]*
      z.design.pairwise.interaction[,col2]
    z.design.names.pairwise.interaction <- c(z.design.names.pairwise.interaction,
                                             #col1-1 is due to there is basline effect in the z design matrix
                                             paste0(tumor.names[col1-1],
                                                    tumor.names[col2-1],
                                                    " interaction effect"))
    z.design.pairwise.interaction <- cbind(z.design.pairwise.interaction,
                                           newcol)
  }
  colnames(z.design.pairwise.interaction) <- z.design.names.pairwise.interaction

  return(z.design.pairwise.interaction)



}


##Generate the z design matrix for saturated model

#' Title
#'
#' @param tumor.character.cat
#' @param tumor.number
#' @param tumor.names
#' @param freq.subtypes
#'
#' @return
#' @export
#'
#' @examples
GenerateZDesignSaturated <- function(tumor.character.cat,tumor.number,tumor.names,freq.subtypes) {
  cutoff <- 10
  z.design.saturated <- GenerateZDesignAdditive(tumor.character.cat,
                                                  tumor.number,
                                                  tumor.names,
                                                  freq.subtypes)
  z.design.names.saturated <- colnames(z.design.saturated)
  ##j represent the order of the interaction
  for(j in 2:tumor.number){
    all.combnation <- combn(tumor.number,j)+1
    combn.number <- ncol(all.combnation)
    #####i represent the ith combination within the jth order interaction
    for(i in 1:combn.number){
      newcol <- rep(1,nrow(z.design.saturated))
      ###k present the kth tumor characteristics within the ith combination
      z.design.names.saturated.one.column <- NULL
      for(k in 1:j){
        col.number <- all.combnation[k,i]
        newcol <-  newcol*z.design.saturated[,col.number]
        z.design.names.saturated.one.column <- paste0(z.design.names.saturated.one.column
                                                      ,tumor.names[col.number-1])
      }
      z.design.names.saturated.one.column <- paste0(z.design.names.saturated.one.column,
                                                    " interaction effect")
      z.design.names.saturated <- c(z.design.names.saturated,
                                    z.design.names.saturated.one.column)
      z.design.saturated <- cbind(z.design.saturated,newcol)

    }

  }
  colnames(z.design.saturated) <- z.design.names.saturated



    return(z.design.saturated)




}
