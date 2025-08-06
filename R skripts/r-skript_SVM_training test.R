###statquest youtube zu SVM


library(limma)
library(statmod)
library(marray)
library(ggplot2)
library(MASS)
library(calibrate)
library(pheatmap)
library(pROC)
library(caret)
library(e1071)
library(caTools)



setwd("/home/tom/Desktop/Analysis/ROC/PDAC")


#Training Set
SVM1 <- read.csv("SVM_input.csv")
head(SVM1)
rownames(SVM1) <- SVM1[,1]
SVM1 <- SVM1[,-1]
head(SVM1)
dim(SVM1)

#Validation/Test set var.1
#set.seed(1234)
split<- sample.split(SVM1$outcome, SplitRatio = 0.7)
training_set<-subset(SVM1,split==TRUE)
test_set<-subset(SVM1, split==FALSE)



                     
#Training SVM
#svm_model <- svm (outcome ~ ., data=training_set, kernel="radial", cost=1, gamma=0.5)
svm_model <- svm (outcome ~ ., data=training_set ,kernel="radial", cost=1, gamma=0.5)

summary(svm_model)


#Testing SVM
#pred2 <- predict(svm_model, test_set[,-ncol(test_set)])
pred2 <- predict(svm_model, test_set[,-ncol(test_set)])

pred2


roc <- roc(test_set$outcome,pred2)
plot(roc,auc.polygon=T,print.thres=T,print.auc=T)

write.csv(pred2,"decision value005.csv")


#Finding best parameters for SVM
#rumprobieren verbesser t evtl ROC kurve
svm_tune <- tune(svm, outcome ~ ., data=training_set, kernel="radial", ranges=list(cost=10^(-3:2), gamma=c(0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05)))
print(svm_tune)

svm_model_tuned <- svm(outcome ~ ., data=training_set, kernel="radial", cost=1, gamma=0.02)

summary(svm_model_tuned)

pred <- predict(svm_model_tuned,test_set[,-ncol(test_set)])
pred

roc <- roc(test_set$outcome,pred)
plot(roc,auc.polygon=T,print.thres=T,print.auc=T) 





