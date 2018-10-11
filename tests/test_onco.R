library(devtools)
library(bc2)
library(readr)





y_standard <- as.matrix(read_csv("./data/onco_y_standard.csv"))
x <- as.matrix(read_csv("./data/gene_and_pc.csv"))
y_onco <- as.matrix(read_csv("./data//onco_y.csv"))
pc <- x[,2:11]
gene <- x[,1]


##############Two stage model MLE
Heter_result = EMmvpoly(y_onco,baselineonly = NULL,main.effect = as.matrix(cbind(x)),pairwise.interaction = NULL,saturated = NULL,missingTumorIndicator = 888)


logodds <- Heter_result[[1]][24:28]
info <- Heter_result[[2]][24:28,24:28]
GlobalTestForAssoc(logodds,info)
GlobalTestForHeter(logodds,info)
IndividualHeterTest(logodds,info)
