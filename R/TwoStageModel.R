

#' Title
#'
#' @param y
#' @param baselineonly
#' @param additive
#' @param pairwise.interaction
#' @param saturated
#' @param missingTumorIndicator
#' @param delta0
#'
#' @return
#' @export
#'
#' @examples
TwoStageModel <- function(y,
                          baselineonly=NULL,
                          additive=NULL,
                          pairwise.interaction=NULL,
                          saturated=NULL,
                          missingTumorIndicator = NULL,
                          delta0 = NULL){
  if(is.null(missingTumorIndicator)==1){
    return(Mvpoly(y,
                  baselineonly,
                  additive,
                  pairwise.interaction,
                  saturated,
                  delta0 = delta0))
  }else{

      return(EMmvpoly(y,
                      baselineonly,
                      additive,
                      pairwise.interaction,
                      saturated,
                      missingTumorIndicator,
                      delta0 = delta0))


  }
}
