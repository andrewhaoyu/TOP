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

