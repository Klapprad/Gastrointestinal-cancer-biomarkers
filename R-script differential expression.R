
################################################### import data ###############################################################################

### set working directory
setwd("/home/tom/Desktop/Analysis/GPR4") #absolute filepath for working directory 
getwd()


### load limma and pheatmap packages
library(limma)
library(pheatmap)
library(pROC)
library(e1071)


### reading and preprocessing microarray data including top and bottom respectively
top <- readTargets("targets_top.txt") #variable top wird definiert
bottom <- readTargets("targets_bottom.txt") #variable bottom wird definiert
dim(top) #dim 
dim(bottom)
head(top)




### data preprocess step1: slides spots quality weight, remove spots which flag is less than or equal to 50 
f <- function(x) as.numeric(x$Flags > -50)
RGtop <- read.maimages(top, source='genepix', wt.fun=f) 
RGbottom <- read.maimages(bottom, source='genepix', wt.fun=f)

### data preprocess step3: merge the data from the top and bottom of one slide



###################################################  data preprocessing  ###########################################################################

### data preprocess step:fast slides quality assessment
# plotMA3by2(RGwhole)


### data preprocess step: background correction methods (normexp or subtract)
RGb.top <- backgroundCorrect(RGtop,method="subtract")
RGb.bottom <- backgroundCorrect(RGbottom,method="subtract")

#RGb <- backgroundCorrect(RGall,method="normexp",offset=50,normexp.method="mle")
#plotDensities(RGb)


### data preprocess step: Within-Array Normalization
MA.top <- normalizeWithinArrays(RGb.top)
MA.bottom <- normalizeWithinArrays(RGb.bottom)
#MA <- normalizeWithinArrays(RGb,method="loess",iterations=10)

plotDensities(MA)

### data preprocess step: Between-Array Normalization
MA.t <- normalizeBetweenArrays(MA.top)
MA.b <- normalizeBetweenArrays(MA.bottom)
plotDensities(MA)
plotMDS(MA)

MA <- rbind(MA.t,MA.b)

### data preprocess step: remove control probes and missing value in gene names
MA <- MA[-grep("control",MA$genes$Name),]
MA <- MA[-which(MA$genes$Name == ""),]
MA <- MA[-which(MA$genes$Name == "empty"),]
MA


### data preprocess step: calculate protein expression value on all the slides by calculating the mean or max of duplicate log2 ratio expression value in one slide
eset <- MA$M
rownames(eset) <- MA$genes$ID
head(eset)
dim(eset)

eset1 <- tapply(eset[,1],INDEX=rownames(eset),FUN=mean,na.rm=T)
head(eset1)
dim(eset1)

eset2 <- as.data.frame(eset1)
colnames(eset2) <- colnames(eset)[1]
head(eset2)
dim(eset2)

for (i in 2:46) {
  eset3 <- tapply(eset[,i],INDEX=rownames(eset),FUN=mean,na.rm=T)
  eset4 <- as.data.frame(eset3)
  colnames(eset4)<- colnames(eset)[i]
  eset2 <- cbind(eset2,eset4)
}
head(eset2)
dim(eset2)

top




design <- modelMatrix(top,ref="Ref")
colnames(design) <- c('CRC_sc','CRC_sr','GC_sc','GC_sr','HCC_sc','HCC_sr','PDAC_sc','PDAC_sr','CTL')
head(design)
dim(design)

##linear model fit
fit <- lmFit(eset2,design)

##construct contrast model
contrast.matrix <- makeContrasts(HCC_sr-PDAC_sr,levels=design)

##differential calculation by contrast model
fit1 <- contrasts.fit(fit,contrast.matrix)

##Bayes test
fit2 <- eBayes(fit1)

##resullt report
topTable(fit2, coef= 1, n=30, adjust="BH")

##plot volcano
volcanoplot(fit2,coef=1, highlight = 6, names=rownames(fit$coefficients), main = "Serum proteins of HCC vs. PDAC")

##result export
dif1 <- topTable(fit2, coef= 1,n = nrow(fit2),adjust="BH")

##screen by p value less than 0.05
dif1p <- dif1[dif1$P.Value <= 0.05,]


dif1pdown <- dif1p[dif1p[,1]<0,]
dif1pup <- dif1p[dif1p[,1]>0,]
dim(dif1pdown)
dim(dif1pup)
head(dif1pup)

#write.csv(dif1, file = "sr_HCC-PDAC.csv")
#write.csv(eset2, file = "mvalues2.csv")


#setwd("/home/tom/Desktop/Analysis/ROC")
# Read in csv files

library(pROC)

roc <- roc(PDAC.csv$outcome,
           PDAC.CSV$BS_089, percent=TRUE,,
           # arguments for auc
           auc=c(100, 90), auc=TRUE,
           auc.focus="sens",
           # arguments for ci
           ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
           # arguments for plot
           plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
           print.auc=TRUE, show.thres=TRUE)