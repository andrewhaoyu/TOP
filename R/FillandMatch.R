
#' Title
#'
#' @param reference
#' @param all
#'
#' @return
#' @export
#'
#' @examples

FillandMatch <- function(reference,all){
  idx.fil <- which(all%in%reference)
  all.fil <- all[idx.fil]
  idx.match <- match(reference,all.fil)
  return(list(idx.fil,idx.match))
}
