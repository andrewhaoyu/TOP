### Generate the frequency of subtypes in the complete data file
#' Title
#'
#' @param y.pheno.complete
#'
#' @return
#' @export
#'
#' @examples
GenerateFreqTable <- function(y.pheno.complete){
  eval_text = NULL
  tumor.number = ncol(y.pheno.complete)-1
  for(i in 2:(tumor.number+1)){
    if(i==(tumor.number+1)){
      eval_text = paste0(eval_text,paste0("y.pheno.complete[,",i,"]"))
    }else{
      eval_text = paste0(eval_text,paste0("y.pheno.complete[,",i,"],"))
    }
  }
  eval_text = paste0("as.data.frame(table(",eval_text,"),stringsAsFactors = F)")
  result = eval(parse(text=eval_text))
  result <- apply(result,2,as.numeric)
  # idx <- which(result[,ncol(result)]==0)
  # result <- result[-idx,,drop=F]
  return(result)
}
