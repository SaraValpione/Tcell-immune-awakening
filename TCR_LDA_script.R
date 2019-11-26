rm(list=ls())

library(rms)
library(mudata)

#You calculate Gini coefficient and Renyi index using LymphoSeq (clonality command) and the script written by Spreafico [Spreafico, R. et al. A circulating reservoir of pathogenic-like CD4+ T cells shares a genetic and phenotypic signature with the inflamed synovial micro-environment. Ann Rheum Dis 75, 459-465, doi:10.1136/annrheumdis-2014-206226 (2016)], then you calculate the delta (W3-T0) for the figure and the quared delta (I called these deltas for ease)for the LDA.
#this below is the analyses for the training and validation cohort, you can download these from GitHub.
#LDA

setwd("/...")
training_deltas<-read.csv("TCR_training_deltas.csv", header=T, sep=",")
validation_deltas<-read.csv("TCR_validation_deltas.csv", header=T, sep=",")

#set the response as a factor
training_deltas$RR<-as.factor(training_deltas$RR)
validation_deltas$RR<-as.factor(validation_deltas$RR)

library(MASS)
fit<-lda(RR~Renyi_PBMC+Gini_PBMC, data=training_deltas,na.action="na.omit")
fit # shows results
plda<-predict(fit, training_deltas)
plda_validation<-predict(fit, validation_deltas)
K<-predict(fit, validation_deltas)$class

library(caret)
modelFit<-train(RR~Renyi_PBMC+Gini_PBMC, method="lda", preProcess=c("scale","center"), data= training_deltas)

prdval<-predict(modelFit, validation_deltas)
confusionMatrix(validation_deltas$RR, predict(modelFit, validation_deltas))

