##############load in the breast cancer example
data(data, package="TOP")
##############this is a simulated breast cancer example
##############there are around 5000 breast cancer cases and 5000 controls, i.e. people without disease
data[1:5,]
##############four different tumor characteristics were included, ER (positive vs negative), PR (positive vs negative), HER2 (positive vs negative), grade (ordinal 1, 2, 3)
##############the phenotype file
y <- data[,1:5]
##############one SNP and one Principal components (PC1) are the covariates
SNP <- data[,6,drop=F]
PC1 <- data[,7,drop=F]
#############fit the additive two-stage model
model.1 <- TwoStageModel(y=y,additive=cbind(SNP,PC1),
                         missingTumorIndicator = 888)
############the model result is a list
############model.1[[4]] are the second stage odds ratio (95% CI) and p-value, the baseline effect is the case-control odds ratio of the reference subtype (ER-,PR-,HER2-,grade0). The main effect are the case-case odds ratio of the tumor characteristics.
model.1[[4]]

############model.1[[5]] are the global association test and global heterogeneity test result of the covariates.
model.1[[5]]
############model.1[[7]] are the case-control odds ratios for all of the subtypes.
model.1[[7]]
############instead of additive model, you can also try different combinations. For example, for the PC1, we use the additive model, but for SNP, we use the baseline only model.
model.2 <- TwoStageModel(y=y,baselineonly = SNP,
                         additive=PC1,
                         missingTumorIndicator = 888)






#load in the breast cancer example
data(data, package="TOP")
#this is a simulated breast cancer example
#there are around 5000 breast cancer cases and 5000 controls, i.e. people without disease
data[1:5,]
#four different tumor characteristics were included, ER (positive vs negative), PR (positive vs negative), HER2 (positive vs negative), grade (ordinal 1, 2, 3)
#the phenotype file
y <- data[,1:5]
#generate the combinations of all the subtypes
#by default, we remove all the subtypes with less than 10 cases
z.standard <- GenerateZstandard(y)
M <- nrow(z.standard) #M is the total number of first stage subtypes

#initial a z.design matrix with M rows, and 5 columns
#each row represent a first stage subtype
#each column represent an aggregated subtype
z.design <- matrix(0,M,5)
#define names for the five intrinsic subtypes
colnames(z.design) <- c("HR+_HER2-_lowgrade",
                        "HR+_HER2+",
                        "HR+_HER2-_highgrade",
                        "HR-_HER2+", 
                        "HR-_HER2-")
#To construct a self design second stage matrix,
#we need to find the correpsonding first stage subtypes
#belonging to specific aggregated subtypes
#for first subtype HR+_HER2-_lowgrade
idx.1 <- which((z.standard[,1]==1|z.standard[,2]==1)
               &z.standard[,3]==0
               &(z.standard[,4]==1|z.standard[,4]==2))
z.design[idx.1,1] <- 1
#for second subtype HR+_HER2+
idx.2 <- which((z.standard[,1]==1|z.standard[,2]==1)
               &z.standard[,3]==1)
z.design[idx.2,2] <- 1
#for third subtype HR+_HER2-_highgrade
idx.3 <- which((z.standard[,1]==1|z.standard[,2]==1)
               &z.standard[,3]==0
               &z.standard[,4]==3)
z.design[idx.3,3] <- 1
#for third subtype HR-_HER2+
idx.4 <- which(z.standard[,1]==0&z.standard[,2]==0
               &z.standard[,3]==1)
z.design[idx.4,4] <- 1
#for third subtype HR-_HER2-
idx.5 <- which(z.standard[,1]==0&z.standard[,2]==0
               &z.standard[,3]==0)
z.design[idx.5,5] <- 1
#one SNP and one Principal components (PC1) are the covariates
SNP <- data[,6,drop=F]
PC1 <- data[,7,drop=F]

model.3 <- EMmvpolySelfDesign(y,
          x.self.design = SNP,
    z.design = z.deign,
    additive=PC1,
  missingTumorIndicator = 888)









############fit the additive two-stage model under the null hypothesis that the second stage parameters of SNP is 0
score.support <- ScoreTestSupportMixedModel(y=y,
                                       additive=PC1,
                                  missingTumorIndicator=888)
  


  









