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
    z.design = z.design,
    additive=PC1,
  missingTumorIndicator = 888)







#fit the two-stage model under the null hypothesis
#that the second stage parameters of SNP is 0
#the model only has one covariate PC1
score.support.fixed <- ScoreTestSupportMixedModel(y=y,
                                                  additive=PC1,
                                                  missingTumorIndicator=888)


#Generate the additive second stage design matrix
z.additive <- cbind(1,z.standard)
#compute the score and information matrix for SNP
score.test.fixed <- ScoreTestMixedModel(y=y,
                                        x=SNP,
                                        z.design = z.additive,
                                        score.test.support=score.support.fixed,
                                        missingTumorIndicator=888)
#the first element is the score
#the second element is the information matrix
score.fixed <- score.test.fixed[[1]]
infor.fixed <- score.test.fixed[[2]]
#compute the global association test p value
DisplayFixedScoreTestResult(score.fixed,infor.fixed) 


  


#we are going to build a two-stage model with 
#baseline parameter and ER case-case parameter as fixed
#We assume the PR, HER2, grade case-case parameter 
#to be random
#fit the two-stage model under the null hypothesis
#that the second stage parameters of SNP is 0
#the model only has one covariate PC1
#Generate the z design matrix for fixed effect
z.design.fixed <- cbind(1,z.standard[,1])
#compute the score and information matrix for fixed effect
score.test.fixed <- ScoreTestMixedModel(y=y,
                    x=SNP,
                    z.design=z.design.fixed,
                    score.test.support=score.support.fixed,
                    missingTumorIndicator=888)
#the first element is the score
#the second element is the information matrix
score.fixed <- score.test.fixed[[1]]
infor.fixed <- score.test.fixed[[2]]

#fit the two-stage model under the null hypothesis
#that only the random effect is 0
#the model will have two covariates, 
#PC1 and the fixed effect of SNP
score.support.random <- ScoreTestSupportMixedModelSelfDesign(y=y,
                        x.self.design  = SNP,
                        z.design = z.design.fixed,
                        additive = PC1,
                        missingTumorIndicator = 888)
#Generate the z design matrix for random effect
#PR, HER2 and grade is random effect
z.design.random <- z.standard[,2:4]

#compute the score and information matrix for random effect
score.test.random <- ScoreTestMixedModel(y=y,
                                         x=SNP,
                                         z.design=z.design.random,
                                         score.test.support=score.support.random,
                                         missingTumorIndicator=888)
#the first element is the score
#the second element is the information matrix
score.random <- score.test.random[[1]]
infor.random <- score.test.random[[2]]

#after we get the the fixed effect score, infor and random effect score, infor, we can combine them through the following function. 
#two p value will be generated.
#the first p value for global association test.
#the second p value is for the null hypothesis 
#that random effect is 0
#Under this situation, the second p value is NOT 
#the global heterogeneity test p value since ER is fixed
p.value.mtop <- DisplayMixedScoreTestResult(
  score.fixed,
  infor.fixed,
  score.random,
  infor.random
)  




