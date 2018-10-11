###Generate the complete y pheno dataframe based on the incomplete data file










#' Title
#'
#' @param y.pheno
#' @param missingTumorIndicator
#'
#' @return
#' @export
#'
#' @examples
GenerateMissingPosition <- function(y.pheno,missingTumorIndicator){


  tumor.number  <-  ncol(y.pheno)-1
  find.missing.position.text = "idx <- which("
  for(i in 2:(tumor.number+1)){
    if(i == (tumor.number+1)){
      find.missing.position.text <- paste0(find.missing.position.text,"y.pheno[,",i,"]==missingTumorIndicator)")
    }else{
      find.missing.position.text <- paste0(find.missing.position.text,"y.pheno[,",i,"]==missingTumorIndicator|")
    }
  }
  eval(parse(text=find.missing.position.text))

  return(idx)




}









#' Title
#'
#' @param y.pheno
#' @param missingTumorIndicator
#'
#' @return
#' @export
#'
#' @examples
GenerateCompleteYPheno <- function(y.pheno,missingTumorIndicator){
  tumor.number  <-  ncol(y.pheno)-1
  find.missing.position.text = "idx <- which("
  for(i in 2:(tumor.number+1)){
    if(i == (tumor.number+1)){
      find.missing.position.text <- paste0(find.missing.position.text,"y.pheno[,",i,"]==missingTumorIndicator)")
    }else{
      find.missing.position.text <- paste0(find.missing.position.text,"y.pheno[,",i,"]==missingTumorIndicator|")
    }
  }
  eval(parse(text=find.missing.position.text))
  if(length(idx)!=0){
    y.pheno.complete = y.pheno[-idx,]
  }else{
    y.pheno.complete = y.pheno
  }

  return(y.pheno.complete)
}



#' Title
#'
#' @param y.pheno
#' @param x.all
#' @param missingTumorIndicator
#'
#' @return
#' @export
#'
#' @examples
GenerateCompleteXCovariates <- function(y.pheno,x.all,missingTumorIndicator){
  tumor.number  <-  ncol(y.pheno)-1
  find.missing.position.text = "idx <- which("
  for(i in 2:(tumor.number+1)){
    if(i == (tumor.number+1)){
      find.missing.position.text <- paste0(find.missing.position.text,"y.pheno[,",i,"]==missingTumorIndicator)")
    }else{
      find.missing.position.text <- paste0(find.missing.position.text,"y.pheno[,",i,"]==missingTumorIndicator|")
    }
  }
  eval(parse(text=find.missing.position.text))
  if(length(idx)!=0){
    x.all.complete <- x.all[-idx,]
  }else{
    x.all.complete <- x.all
  }
  return(x.all.complete)
}

