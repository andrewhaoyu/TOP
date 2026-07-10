###Generate the complete y pheno dataframe based on the incomplete data file










#' Title
#'
#' @param y.pheno Input for `y.pheno`.
#' @param missingTumorIndicator Input for `missingTumorIndicator`.
#'
#' @return A result from this computation.
#' @export
#'
GenerateMissingPosition <- function(y.pheno,missingTumorIndicator){
  tumor <- y.pheno[, -1, drop = FALSE]
  matches <- !is.na(tumor) & tumor == missingTumorIndicator
  which(rowSums(matches) > 0L)




}









#' Title
#'
#' @param y.pheno Input for `y.pheno`.
#' @param missingTumorIndicator Input for `missingTumorIndicator`.
#'
#' @return A result from this computation.
#' @export
#'
GenerateCompleteYPheno <- function(y.pheno,missingTumorIndicator){
  idx <- GenerateMissingPosition(y.pheno, missingTumorIndicator)
  if (length(idx)) y.pheno[-idx, , drop = FALSE] else y.pheno
}



#' Title
#'
#' @param y.pheno Input for `y.pheno`.
#' @param x.all Input for `x.all`.
#' @param missingTumorIndicator Input for `missingTumorIndicator`.
#'
#' @return A result from this computation.
#' @export
#'
GenerateCompleteXCovariates <- function(y.pheno,x.all,missingTumorIndicator){
  idx <- GenerateMissingPosition(y.pheno, missingTumorIndicator)
  if (length(idx)) x.all[-idx, , drop = FALSE] else x.all
}
