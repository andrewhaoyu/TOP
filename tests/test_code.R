library(TOP)
#load in the breast cancer example
data(data, package="TOP")
#this is a simulated breast cancer example
#there are around 5000 breast cancer cases and 5000 controls disease
data[1:5,]
#four different tumor characteristics were included, 
#ER (positive vs negative), 
#PR (positive vs negative),
#HER2 (positive vs negative)
#grade (oridinal 1,2,3)
#the phenotype file
y <- data[,1:5]
#generate the combinations of all the subtypes
#by default, we remove all the subtypes with less than 10 cases
z.standard <- GenerateZstandard(y)
M <- nrow(z.standard) #M is the total number of first stage subtypes
y=y;
additive=PC1;
missingTumorIndicator=888;
baselineonly=NULL;
pairwise.interaction=NULL;
saturated=NULL;
delta0 = NULL;
cutoff =10


data(data, package="TOP")
#this is a simulated breast cancer example
#there are around 5000 breast cancer cases and 5000 controls disease
data[1:5,]
#four different tumor characteristics were included, 
#ER (positive vs negative), 
#PR (positive vs negative),
#HER2 (positive vs negative)
#the phenotype file
y <- data[,1:5]
#one SNP
#one Principal components (PC1) are the covariates
SNP <- data[,6,drop=F]
PC1 <- data[,7,drop=F]
#fit the additive two-stage model
model.1 <- TwoStageModel(y=y,additive=cbind(SNP,PC1),
                         missingTumorIndicator = 888,
                         cutoff=10)
temp1 = model.1[[1]]

model.1 <- TwoStageModel(y=y.pheno.complete,additive=cbind(SNP,PC1),
                         missingTumorIndicator = NULL,
                         cutoff=10)
temp2 = model.1[[1]]




y=y.c;
baselineonly=NULL;
additive=PC1.c;
pairwise.interaction=NULL;
saturated=NULL;
missingTumorIndicator = NULL;
delta0= NULL;
cutoff = 10;


ScoreTestSupportMixedModel <- function(y,
                                       baselineonly=NULL,
                                       additive=NULL,
                                       pairwise.interaction=NULL,
                                       saturated=NULL,
                                       missingTumorIndicator = 888,
                                       delta0 = NULL,
                                       cutoff =10){
  
  y <- as.matrix(y)
  tumor.number <- ncol(y)-1
  y.case.control <- y[,1,drop=F]
  y.tumor <- y[,2:(tumor.number+1),drop=F]
  if(is.null(missingTumorIndicator)==T){
    y.pheno.complete = y  
  }else{
    y.pheno.complete <- GenerateCompleteYPheno(y,missingTumorIndicator)
  }
  
  freq.subtypes <- GenerateFreqTable(y.pheno.complete)
  if(CheckControlTumor(y.case.control,y.tumor)==1){
    return(print("ERROR:The tumor characteristics for control subtypes should put as NA"))
  }
  tumor.names <- colnames(y.tumor)
  if(is.null(tumor.names)){
    tumor.names <- paste0(c(1:tumor.number))
  }
  tumor.character.cat = GenerateTumorCharacterCat(y.pheno.complete)
  z.design.baselineonly <- GenerateZDesignBaselineonly(tumor.character.cat,
                                                       tumor.number,
                                                       tumor.names,
                                                       freq.subtypes,
                                                       cutoff)
  z.design.additive <- GenerateZDesignAdditive(tumor.character.cat,
                                               tumor.number,
                                               tumor.names,
                                               freq.subtypes,
                                               cutoff)
  
  if(tumor.number>=2){
    z.design.pairwise.interaction <- GenerateZDesignPairwiseInteraction(tumor.character.cat,
                                                                        tumor.number,
                                                                        tumor.names,
                                                                        freq.subtypes,
                                                                        cutoff)
    z.design.saturated <- GenerateZDesignSaturated(tumor.character.cat,
                                                   tumor.number,
                                                   tumor.names,
                                                   freq.subtypes,
                                                   cutoff)
    
  }else{
    z.design.pairwise.interaction <- z.design.additive
    z.design.saturated <- z.design.additive
    
  }
  z.all <- ZDesigntoZall(baselineonly,
                         additive,
                         pairwise.interaction,
                         saturated,
                         z.design.baselineonly,
                         z.design.additive,
                         z.design.pairwise.interaction,
                         z.design.saturated)
  
  if(is.null(delta0)==T){
    delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
  }else{
    delta0 =delta0
  }
  
  #x.all has no intercept yet
  #we will add the intercept in C code
  x.all <- GenerateXAll(y,baselineonly,additive,pairwise.interaction,saturated)
  ###z standard matrix means the additive model z design matrix without baseline effect
  ###z standard matrix is used to match the missing tumor characteristics to the complete subtypes
  
  z.standard <- z.design.additive[,-1,drop=F]
  
  Score.Support = EMStepScoreTestSupportMixedModel(delta0,y,x.all,z.standard,z.all,missingTumorIndicator)
  
  
  # score_support_result <- score_support(pxx,x.all,baselineonly,z.all,z.standard,y_em)
  #score_test_mis <- score_test_mis(y_em,baselineonly,score_support_result)
  #return(list(score_c=score_test_mis$score_c,infor_c = score_test_mis$infor_c))
  result <- Score.Support
  result[[7]] <- z.standard
  return(result)
  
}



