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
x.self.design= SNP.c;
z.design = z.design;
additive=PC1.c;
pairwise.interaction=NULL;
saturated=NULL;
missingTumorIndicator = NULL;
delta0= NULL;
cutoff = 10;


y.pheno.complete <- y


freq.subtypes <- GenerateFreqTable(y.pheno.complete)
tumor.number <- ncol(y)-1
y.case.control <- y[,1,drop=F]
y.tumor <- y[,2:(tumor.number+1),drop=F]
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
                                             freq.subtypes)
if(tumor.number>=2){
  z.design.pairwise.interaction <- GenerateZDesignPairwiseInteraction(tumor.character.cat,
                                                                      tumor.number,
                                                                      tumor.names,
                                                                      freq.subtypes)
  z.design.saturated <- GenerateZDesignSaturated(tumor.character.cat,
                                                 tumor.number,
                                                 tumor.names,
                                                 freq.subtypes)
  
}else{
  z.design.pairwise.interaction <- z.design.additive
  z.design.saturated <- z.design.additive
  
}
full.second.stage.names <- colnames(z.design.saturated)
covar.names <- GenerateCovarName(baselineonly,
                                 additive,
                                 pairwise.interaction,
                                 saturated,
                                 x.self.design = x.self.design)
z.all <- ZSelfDesigntoZall(x.self.design,
                           baselineonly,
                           additive,
                           pairwise.interaction,
                           saturated,
                           z.design,
                           z.design.baselineonly,
                           z.design.additive,
                           z.design.pairwise.interaction,
                           z.design.saturated)
delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
z.standard <- z.design.additive[,-1,drop=F]
if(is.null(delta0)){
  delta0 <-StartValueFunction(freq.subtypes,y.case.control,z.all)
}

#x.all has no intercept yet
#we will add the intercept in C code
x.all <- GenerateSelfXAll(y,x.self.design,baselineonly,additive,pairwise.interaction,saturated)
z.standard <- z.design.additive[,-1,drop=F]


fit.result <- ProbFitting(delta0,y,x.all,z.standard,z.all,missingTumorIndicator=888)
y.fit <- fit.result[[1]]
idx.drop = fit.result[[4]]

if(length(idx.drop)!=0){
  x.all <- x.all[-idx.drop,,drop=F]
  y.fit <- y.fit[-idx.drop,,drop=F]
}




x.all <- as.matrix(x.all)
z.standard <- z.design.additive[,-1]
M <- as.integer(nrow(z.standard))

tol <- as.numeric(1e-04)


delta_old <- delta0

N <- as.integer(nrow(x.all))

M <- as.integer(nrow(z.standard))

NCOV   <- as.integer(ncol(x.all))
NM     <- N*M
nparm  <- as.integer(length(delta0))
deltai <- as.numeric(delta0)

NITER  <- as.integer(500)
Y <- as.numeric(as.vector(y.fit))
X <- as.numeric(as.vector(x.all))
ZallVec = as.numeric(as.vector(z.all))
Znr = as.integer(nrow(z.all))
Znc = as.integer(ncol(z.all))
debug     <- as.integer(1)
ret_rc    <- as.integer(1)
ret_delta <- as.numeric(rep(-9999, nparm))
ret_info <- as.numeric(rep(-9999,nparm^2))
ret_p <- as.numeric(rep(0,NM))
ret_lxx <- as.numeric(rep(0,NM))
loglikelihood <- as.numeric(-1);



temp <- .C("Mvpoly_complete",deltai, nparm, Y=Y, X, ZallVec,Znr,Znc, N, M, NCOV, NITER, tol,
           debug, ret_rc=ret_rc, ret_delta=ret_delta,ret_info=ret_info,ret_p=ret_p,loglikelihood = loglikelihood)

info <- matrix(unlist(temp$ret_info),nparm,nparm)
result <- list(temp$ret_delta,info,
               temp$ret_p)


# infor_mis_c <- infor_mis(y_em,x.all,z.all)
#infor_obs <- result[[2]]-infor_mis_c
delta=result[[1]]
infor_obs=result[[2]]
p=result[[3]]
loglikelihood = temp$loglikelihood
AIC = 2*nparm - 2*loglikelihood


Mvpoly.result <- (list(delta=delta,
                       infor_obs=infor_obs,
                       p=p,y_em=NULL,
                       M=M,
                       NumberofTumor=ncol(z.standard),
                       loglikelihood = loglikelihood,
                       AIC = AIC
))
###delta represent second stage parameters
delta <- Mvpoly.result$delta
covariance.delta <- solve(Mvpoly.result$infor_obs)
loglikelihood <- Mvpoly.result$loglikelihood
AIC <- Mvpoly.result$AIC
second.stage.mat <-
  GenerateSecondStageMat(baselineonly,
                         additive,
                         pairwise.interaction,
                         saturated,
                         M,
                         full.second.stage.names,
                         covar.names,
                         delta,
                         z.design.additive,
                         z.design.pairwise.interaction,
                         z.design.saturated)
##take out the intercept from second stage parameters

takeout.intercept.result <- TakeoutIntercept(delta,covariance.delta,
                                             M,
                                             tumor.names,
                                             z.all,covar.names)
beta <- takeout.intercept.result$beta
covariance.beta <- takeout.intercept.result$covariance.beta
delta.no.inter <- takeout.intercept.result$delta.no.inter
covariance.delta.no.inter <-
  takeout.intercept.result$covariance.delta.no.inter
beta.no.inter <- takeout.intercept.result$beta.no.inter
covariance.beta.no.inter <- takeout.intercept.result$covariance.beta.no.inter



second.stage.test <- SecondStageTest(delta.no.inter,covariance.delta.no.inter,M,second.stage.mat)
global.test <- GenerateGlobalTest(delta.no.inter,
                                  covariance.delta.no.inter,
                                  M,
                                  second.stage.mat)
##beta represent first stage parameters

subtypes.names <- GenerateSubtypesName(z.design.additive,M,
                                       tumor.names)
first.stage.mat <- GenerateFirstStageMat(beta,
                                         covar.names,
                                         subtypes.names)

first.stage.test <- FirstStageTest(beta.no.inter,
                                   covariance.beta.no.inter,
                                   M,
                                   first.stage.mat)



return(list(delta=delta,covariance.delta=covariance.delta,second.stage.mat = second.stage.mat,second.stage.test,global.test,first.stage.mat,first.stage.test,loglikelihood = loglikelihood,
            AIC = AIC,beta=beta,covariance.beta=covariance.beta,
            z.standard=z.standard))