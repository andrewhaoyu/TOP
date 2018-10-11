#' Title
#'
#' @param freq.subtypes
#' @param y.case.control
#' @param z.all
#'
#' @return
#' @export
#'
#' @examples
StartValueFunction = function(freq.subtypes,y.case.control,z.all){
  ###cutoff for take one subject
  cutoff=10
  ncontrol <- sum(y.case.control==0)
  p.freq <- ncol(freq.subtypes)
  freq = freq.subtypes[,p.freq]
  idx =which(freq<=cutoff)
  if(length(idx)!=0){
    freq.subtypes = freq.subtypes[-idx,]
    freq = freq.subtypes[,p.freq]
    total = sum(freq)+ncontrol
    p.empirical = freq/total
    delta_inter = log(p.empirical/(1-p.empirical))
  }else{
    total = sum(freq)+ncontrol
    p.empirical = freq/total
    delta_inter = log(p.empirical/(1-p.empirical))
  }
  delta0 <- rep(0,ncol(z.all))
  delta0[1:length(delta_inter)] <- delta_inter
  return(delta0)
}
