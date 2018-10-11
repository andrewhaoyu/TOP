
#' Title
#'
#' @param logodds
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
GlobalTestForAssoc <- function(logodds,sigma){
  sigma <- as.matrix(sigma)
  df <- length(logodds)
  GTA.stat <- t(logodds)%*%solve(sigma)%*%logodds
  p.value.GTA <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)
  places = 3
  power.number <- floor(-log10(p.value.GTA))+places
  ###format the output with three digits in total
  p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)

  return(p.value.GTA)

}





#' Title
#'
#' @param logodds
#' @param sigma
#' @param self.design
#'
#' @return
#' @export
#'
#' @examples
GlobalTestForHeter <- function(logodds,sigma,self.design=F){

  if(self.design==F){
    if(length(logodds)==1){
      return(NA)
    }else{
      sigma <- as.matrix(sigma)
      df <- length(logodds)
      sigma.casecase <- sigma[2:df,2:df]
      logodds.casecase <- logodds[2:df]
      GTH.stat <- t(logodds.casecase)%*%solve(sigma.casecase)%*%logodds.casecase
      p.value.GTH <- pchisq(as.numeric(GTH.stat),df=(df-1),lower.tail = F)
      places <- 3
      power.number <- floor(-log10(p.value.GTH))+places
      p.value.GTH <- round(p.value.GTH*10^power.number)/(10^power.number)


      return(p.value.GTH)
    }
    }else{
      if(length(logodds)==1){
        return(NA)
      }else{

        sigma <- as.matrix(sigma)
        df <- length(logodds)
        z.design <- diag(df-1)
        z.design <- cbind(-1,z.design)


        logodds.casecase <- z.design%*%logodds
        sigma.casecase <- z.design%*%sigma%*%t(z.design)
        GTH.stat <- t(logodds.casecase)%*%solve(sigma.casecase)%*%logodds.casecase
        p.value.GTH <- pchisq(as.numeric(GTH.stat),df=(df-1),lower.tail = F)
        places <- 3
        power.number <- floor(-log10(p.value.GTH))+places
        p.value.GTH <- round(p.value.GTH*10^power.number)/(10^power.number)


        return(p.value.GTH)
    }

  }

}


#' Title
#'
#' @param logodds
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
IndividualHeterTest <- function(logodds,sigma){

  sigma <- as.matrix(sigma)
  var.logodds <- diag(sigma)
  df <- length(logodds)
  z <- logodds/sqrt(var.logodds)
  p.value.IHT <- PvalueFunction(z)

  return(p.value.IHT)

}


#' Title
#'
#' @param score
#' @param infor
#'
#' @return
#' @export
#'
#' @examples
ScoreGlobalTestForAssoc <- function(score,infor){
  infor <- as.matrix(infor)
  df <- length(score)
  GTA.stat <- score%*%solve(infor)%*%t(score)
  p.value.GTA <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)
  places = 3
  power.number <- floor(-log10(p.value.GTA))+places
  ###format the output with three digits in total
  p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)

  return(p.value.GTA)

}


#' Title
#'
#' @param score
#' @param infor
#'
#' @return
#' @export
#'
#' @examples
ScoreGlobalMixedTestForAssoc <- function(score,infor){
  infor <- as.matrix(infor)
  df <- length(score)
  GTA.stat <- score%*%solve(infor)%*%t(score)
  p.value.GTA <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)
  places = 3
  power.number <- floor(-log10(p.value.GTA))+places
  ###format the output with three digits in total
  p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)

  GTA.mixed <- score%*%t(score)

  lamda <- eigen(infor)$values

  result <- davies(GTA.mixed,lamda,lim = 2000000,acc=1e-9)
  p.value.GTA.mixed <- result[[3]]

  if(result[[2]]!=0){
    print("chisq p value accuracy could't be reached")
  }

  if(p.value.GTA.mixed <0){
    p.value.GTA.mixed <- 1e-09
  }



  return(c(p.value.GTA,p.value.GTA.mixed))

}


#' Title
#'
#' @param score.casecase
#' @param infor.casecase
#'
#' @return
#' @export
#'
#' @examples
ScoreMixedGlobalTestForHeter <- function(score.casecase,infor.casecase){


    GTH.stat <- as.numeric(score.casecase%*%t(score.casecase))
    lamda <- eigen(infor.casecase)$values

    result <- davies(GTH.stat,lamda,lim = 2000000,acc=1e-9)
    p.value.GTH <- result[[3]]

    if(result[[2]]!=0){
      print("chisq p value accuracy could't be reached")
    }

    if(p.value.GTH <0){
      p.value.GTH <- 1e-09
    }
    #places <- 3
    #power.number <- floor(-log10(p.value.GTH))+places
    #p.value.GTH <- round(p.value.GTH*10^power.number)/(10^power.number)


    return(p.value.GTH)



}

#' Title
#'
#' @param p.value.score.heter
#' @param score.baseline
#' @param infor.baseline
#'
#' @return
#' @import CompQuadForm
#' @export
#'
#' @examples
ScoreMixedGlobalTestForAssoc <- function(p.value.score.heter,
                                         score.baseline,
                                         infor.baseline){
  infor.baseline <- as.matrix(infor.baseline)
  df <- length(score.baseline)
  GTA.stat <- score.baseline%*%solve(infor.baseline)%*%t(score.baseline)
  p.value.baseline <- pchisq(as.numeric(GTA.stat),df=df,lower.tail = F)

  mix.stat <- -2*log(p.value.baseline)-2*log(p.value.score.heter)

  p.value.mixed <- pchisq(as.numeric(mix.stat),df=4,lower.tail = F)

  return(p.value.mixed)




  #places = 3
  #power.number <- floor(-log10(p.value.GTA))+places
  ###format the output with three digits in total
  #p.value.GTA <- round(p.value.GTA*10^power.number)/(10^power.number)

}

