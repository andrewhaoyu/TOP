#' Generate the three Z design matrix given the tumor characteristics
#'
#' @param y the phenotype file. The first column is the case control disease status. The other columns are the tumor characteristics status
#'
#' @return
#' @export
#'
#' @examples
GenerateZDesignCombination <- function(y,
                                       missingTumorIndicator = 888,
                                       cutoff = 10){
  
  missing.data.vec <- GenerateMissingPosition(y,missingTumorIndicator)
  y.pheno.complete <- y[-missing.data.vec,,drop=F]
  initial.set <- InitialSetup(y.pheno.complete,
                              baselineonly=NULL,
                              additive=NULL,
                              pairwise.interaction=NULL,
                              saturated=NULL,
                              cutoff =cutoff
  )
  
  
  z.design.additive = initial.set[[5]]
  z.design.pairwise.interaction = initial.set[[6]]
  z.design.saturated = initial.set[[7]]
  return(list(z.design.additive = z.design.additive,
              z.design.pairwise.interaction = z.design.pairwise.interaction,
              z.design.saturated = z.design.saturated
  ))
}
