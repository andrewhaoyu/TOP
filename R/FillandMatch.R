
#' Title
#'
#' @param reference Input for `reference`.
#' @param all Input for `all`.
#'
#' @return A result from this computation.
#' @export
#'

FillandMatch <- function(reference,all){
  idx.fil <- which(all%in%reference)
  all.fil <- all[idx.fil]
  idx.match <- match(reference,all.fil)
  return(list(idx.fil,idx.match))
}
