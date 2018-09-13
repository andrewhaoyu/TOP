#' Title
#'
#' @param z
#'
#' @return
#' @export
#'
#' @examples
PvalueFunction <- function(z){
  result <- NULL
  for(i in 1:length(z)){

    result <- c(result,2*(pnorm(-abs(z[i]))))

  }
  places = 3
  ####round p.value under scientific accuracy
  accuracy <- floor(-log10(result))+places
  result <- round(result*10^accuracy)/(10^accuracy)


  return(result)
}
