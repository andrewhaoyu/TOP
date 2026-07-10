

#' Title
#'
#' @param y Input for `y`.
#' @param baselineonly Input for `baselineonly`.
#' @param additive Input for `additive`.
#' @param pairwise.interaction Input for `pairwise.interaction`.
#' @param saturated Input for `saturated`.
#' @param missingTumorIndicator Input for `missingTumorIndicator`.
#' @param cutoff Subtypes with counts at or below this value are removed.
#'
#' @return A result from this computation.
#' @export
#'
ScoreTestSupport <- function(y,
                               baselineonly=NULL,
                               additive=NULL,
                               pairwise.interaction=NULL,
                               saturated=NULL,
                               missingTumorIndicator = NULL,
                               cutoff = 10
                             ){
 if(is.null(missingTumorIndicator)==1){
  return(CompleteCasesScoreTestSupport(y,
                                                   baselineonly,
                                                   additive,
                                                   pairwise.interaction,
                                                   saturated,
                                                   cutoff = cutoff))
 }else{
   #if(missingDataAlgorithm=="EM"){
    return(EMScoreTestSupport(y,
                                          baselineonly,
                                          additive,
                                          pairwise.interaction,
                                          saturated,
                                          missingTumorIndicator,
                                          cutoff = cutoff))
  # }else if(missingDataAlgorithm=="OneStepMLE"){

   #}
 }

}
