### Generate the potential tumor characteristic category(binary or categorical)
#' Title
#'
#' @param y.pheno.complete
#'
#' @return
#' @export
#'
#' @examples
GenerateTumorCharacterCat <- function(y.pheno.complete){
  tumor.number = ncol(y.pheno.complete)-1
  y.tumor.complete = y.pheno.complete[,2:(tumor.number+1)]
  tumor.character.cat = list()
  for(i in 1:tumor.number){
    unique.tumor.cat = unique(y.tumor.complete[!is.na(y.tumor.complete[,i]),i])
    unique.tumor.cat = unique.tumor.cat[order(unique.tumor.cat)]
    tumor.character.cat[[i]] = unique.tumor.cat
  }
  return(tumor.character.cat)
}
