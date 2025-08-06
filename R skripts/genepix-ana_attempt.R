############################## harriet's attempt at genepix analysis ##############################
# i am sorry for what you are about to witness

# before you start ----

# any section with "CHANGE:" above it means you need to change some detail for the code to work (eg setwd)
# "=>" is what I saved the plot / data as
# ">  " shows what is returned when i do it (incl error/warning messages etc)

# setting up ----

# CHANGE: wd
setwd("/home/tom/Desktop/Analysis/GPR2")
getwd()

## downloading libraries ----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("marray")
install.packages("calibrate")

intall.packages("limma")
library(limma)
# needed

install.packages("pheatmap")
library(pheatmap)
# ?

install.packages("pROC")
library(pROC)
# ?

install.packages("e1071")
library(e1071)
# ?

install.packages("statmod")
library(statmod)
# needed

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


# reading microarray data, top and bottom separately ----

# CHANGE: titles of your target files if they aren't titled like mine are
top <- readTargets("targets_top.txt")
bottom <- readTargets("targets_bottom.txt")

# def: wt.fun -> function to calculate spot quality weights
# source = 'genepix' tells the function what the setup of the thing it is reading is
f <- function(x) as.numeric(x$Flags > -50)
RGtop <- read.maimages(top, source='genepix', wt.fun=f) 
RGbottom <- read.maimages(bottom, source='genepix', wt.fun=f)
# note: RG top and bottom are LISTS

# sidenote: if your array slide has a different # of columns / rows do the following:
#RG$printer$ngrid.c <- #COLUMNS 
# make sure you do it in the correct RG file!

# before preprocessing ----

# combining the top and bottom for no purpose other than the occasional diagram bc most functions don't work on the combined file
RGwhole <- rbind(RGtop,RGbottom)
RGwhole$printer$ngrid.r <- 12 # the number of rows of blocks on the slide (6 on top + 6 on bottom)
str(RGwhole) # just shows us that changing the ngrid.r worked

## plotting some things ----
plotDensities(RGwhole) #basic histogram-y thing w the red and the green values of the slides
plotDensities(RGtop)
plotDensities(RGbottom)
# this is what will happen if you imageplot on the RGwhole:
  # imageplot(log2(RGwhole$Rb[,1]), RGwhole$printer, low="white", high="red")
    # > Error: Number of image spots does not agree with layout dimensions
# which is why we do everything separately
# do the $printer thing to determine the # of grid rows and columns on the array
RGtop$printer <- getLayout(RGtop$genes)
imageplot(log2(RGtop$Rb[,1]), RGtop$printer, low = "white", high = "red") # => RGtop_Rb_imageplot in b4_normie
RGbottom$printer <- getLayout(RGbottom$genes)
imageplot(log2(RGbottom$Rb[,1]), RGbottom$printer, low = "white", high = "red") # => RGbottom_Rb_imageplot in b4_normie
# same thing but w Gb (green background)
imageplot(log2(RGtop$Gb[,1]), RGtop$printer, low="white", high="green") # => RGtop_Gb_imageplot in b4_normie
imageplot(log2(RGbottom$Gb[,1]), RGbottom$printer, low="white", high="green") # => RGbottom_Gb_imageplot in b4_normie
# imageplot shows spatial effects across the microarray, like some regions being more fluorescent(?) for some reason
# you can use the argument main="desired_title" to change the the title of the diagram

## boxplots. why? idk. ----
boxplot(data.frame(log2(RGwhole$Gb)),main="Green background", las = 2) # => RGwhole_Gb_boxplot
# > Warning message:
# > In bplt(at[i], wid = width[i], stats = z$stats[, i], out = z$out[z$group ==  :
  # > Outlier (-Inf) in boxplot 2 is not drawn
# ^^^idk what that means but i did get a plot........
boxplot(data.frame(log2(RGwhole$Rb)),main="Red background", las = 2) # => RGwhole_Rb_boxplot

imageplot(log2(RGtop$Gb[,1]),RGtop$printer) # => RGtop_blue_imageplot
imageplot(log2(RGbottom$Gb[,1]),RGbottom$printer) # => RGbottom_blue_imageplot
# ??? idk why we wanted those but i did them!

# the pre-processing ----

# note : the x is what i used to note that the data were corrected

## background correction ----
RGtopx <- backgroundCorrect(RGtop, method = "normexp", offset = 50)
# both fawaz and chaoyang had the offset = 50 sooo i am doing that too
plotDensities(RGtopx) # => RGtopx_densities
# plotDensities is another histogram of expression values
RGbottomx <- backgroundCorrect(RGbottom, method="normexp", offset = 50)
plotDensities(RGbottomx) # => RGbottomx_densities
RGwholex <- backgroundCorrect(RGwhole, method="normexp", offset = 50)
plotDensities(RGwholex) # => RGwholex_densities

## normalize within the arrays ----
# MA = microarray
# this is one of the things where the "whole" dataset won't work
# MAwhole <- normalizeWithinArrays(RGwholex, method = "loess")
# plotPrintTipLoess(MAwhole)
  #  > Error: arguments imply differing number of rows: 5808, 6336
MAtop <- normalizeWithinArrays(RGtopx, method = "loess")
plotPrintTipLoess(MAtop) # => MAtop_PrintTipLoess
# idk what the plot means but there are dots!! and lines!!!
MAbottom <- normalizeWithinArrays(RGbottomx, method = "loess")
plotPrintTipLoess(MAbottom) # => MAbottom_PrintTipLoess
# boxplots too apparently
boxplot(MAtop$M~col(MAtop$M), names=colnames(MAtop$M), cex.lab = 0.75, las = 2, cex.axis = 0.55, main = "MAtop M-Values") # => MAtop_Mvalues_boxplot
# boxplot(quantitative_variable ~ qualitative_variable)
# so it makes a boxplot of the data, and splits the data based on the qualitative variable for side-by-side boxplots
# in this case, M components of MAtop (M is a matrix of the log2 expression values)
# col returning a matrix of integers indicating their column number
# cex.lab, las, cex.axis, and main are just me making the plot pretty
boxplot(MAbottom$M~col(MAbottom$M), names=colnames(MAbottom$M), cex.lab = 0.65, las = 1, cex.axis = 0.5, main = "MAbottom M-Values") # => MAbottom_Mvalues_boxplot

## now we normalize between arrays ----
MAtopx <- normalizeBetweenArrays(MAtop,method="scale")
boxplot(MAtopx$M~col(MAtopx$M), names=colnames(MAtopx$M), cex.lab = 0.75, las = 2, cex.axis = 0.55, main = "MAtop M-Values, norm bw arrays") # => MAtopx_Mvalues_normbw_boxplot
MAbottomx <- normalizeBetweenArrays(MAbottom,method="scale")
boxplot(MAbottomx$M~col(MAbottomx$M), names=colnames(MAbottomx$M), cex.lab = 0.65, las = 1, cex.axis = 0.5, main = "MAbottom M-Values, norm bw arrays") # => MAbottomx_Mvalues_normbw_boxplot
# again, the boxplot doesn't work on the "whole" data

## removing controls, blanks, and empties ----
### controls ----
MAtopx <- MAtopx[-grep("control", MAtopx$genes$Name),]
MAbottomx <- MAbottomx[-grep("control", MAbottomx$genes$Name),]
### blanks ----
MAtopx <- MAtopx[-which(MAtopx$genes$Name == ""),]
MAbottomx <- MAbottomx[-which(MAbottomx$genes$Name == ""),]
### empties ----
MAtopx <- MAtopx[-which(MAtopx$genes$Name == "empty"),]
MAbottomx <- MAbottomx[-which(MAbottomx$genes$Name == "empty"),]

## calculating protein exp value ----
# "by calculating the mean of duplicate log2 ratio expression value in one slide"
# eset is set of expression values
# M is a matrix of the log2 expression values
eset_top <- MAtopx$M
rownames(eset_top) <- MAtopx$genes$ID
head(eset_top)
# all of the gene names are doubled??? should it be like that? i do not know
dim(eset_top)
# > [1] 2384  8
# gene names as rows, samples as columns. these are the expression values. i have 8 slides so that makes sense
# so now we apply a function (mean) over each value in each column
eset_top1 <- tapply(eset_top[,1], INDEX = rownames(eset_top), FUN = mean, na.rm = T)
head(eset_top1)
# so the genes that were in the head of the eset_top[,1] are different to those in the head of eset_top1....
# should it be like that???
# also eset_top1 is a "double" whatever the hell that means
dim(eset_top1)
# > [1] 1179
#repeating all of that with the bottom
eset_bottom <- MAbottomx$M
rownames(eset_bottom) <- MAbottomx$genes$ID
head(eset_bottom)
dim(eset_bottom)
# > [1] 3018  8
eset_bottom1 <- tapply(eset_bottom[,1], INDEX = rownames(eset_bottom), FUN = mean, na.rm = T)
head(eset_bottom1)
# same thing abt how the heads of eset_bottom1 and eset_bottom[,1] differ....
dim(eset_bottom1)
# > [1] 1498
# next turn the eset1 into a dataframe and rename the columns. you are setting up for the upcoming for loop
eset_top2 <- as.data.frame(eset_top1)
colnames(eset_top2) <- colnames(eset_top)[1]
head(eset_top2)
dim(eset_top2)
# > [1] 1179  1
eset_bottom2 <- as.data.frame(eset_bottom1)
colnames(eset_bottom2) <- colnames(eset_bottom)[1]
head(eset_bottom2)
dim(eset_bottom2)
# > [1] 1498  1
# YOU NEED TO HAVE THE ESET_2s PRE-EXISTING BECAUSE OTHERWISE IT OVERWRITES THE VALUES SO YOU ONLY END UP W THE RESULTS FROM THE ONE COLUMN
# hence the cbind in the for-loop
for (i in 2:ncol(eset_top)) {
  eset_top3 <- tapply(eset_top[,i], INDEX = rownames(eset_top), FUN = mean, na.rm = T)
  eset_top4 <- as.data.frame(eset_top3)
  colnames(eset_top4) <- colnames(eset_top)[i]
  eset_top2 <- cbind(eset_top2, eset_top4)
}
dim(eset_top2)
# > [1] 1179   8
# please note the i in 2:ncol etc is TWO BECAUSE YOU ALREADY HAVE THE FIRST COLUMN DONE THAT'S WHAT YOU DID TO SET IT UP
head(eset_top2)
#same thing for the bottom (you already have done eset_bottom2 sooo just continue from the for-loop)
for (i in 2:ncol(eset_bottom)) {
  eset_bottom3 <- tapply(eset_bottom[,i], INDEX = rownames(eset_bottom), FUN = mean, na.rm = T)
  eset_bottom4 <- as.data.frame(eset_bottom3)
  colnames(eset_bottom4) <- colnames(eset_bottom)[i]
  eset_bottom2 <- cbind(eset_bottom2, eset_bottom4)
}
dim(eset_bottom2)
# > [1] 1498   8
head(eset_bottom2)
# all looks well! ...for now...


########### hier daten zusammenführen damit ein Blot
#combined <- rbind(eset_top, eset_bottom)
#combined


## now to save all the results so far. saves them in your WD ----
write.csv(eset_top2, file = "mvalue_top.csv")
write.csv(eset_bottom2, file = "mvalue_bottom.csv")

## get your array weights. not sure how this is useful. ----
arrayw_top <- arrayWeights(MAtopx)
barplot(arrayw_top, xlab="Array", ylab="Weight", col="white", las=2, main = "Array Weights Top")
abline(h=1, lwd=1, lty=2)
# => array_weights_top
arrayw_bottom <- arrayWeights(MAbottomx)
barplot(arrayw_bottom, xlab="Array", ylab="Weight", col="white", las=2, main = "Array Weights Bottom")
abline(h=1, lwd=1, lty=2)
# => array_weights_bottom

# now we do the analysis ----

## generate design matrix: ----
# modelMatrix makes a design matrix from RNA target information for a 2-colour microarray exp (using your targets files!)
# CHANGE : names of your target files! mine are "top" and "bottom"
design_top <- modelMatrix(targets_top, ref="Ref")
# > Found unique target names:
  # > Cancer PSC Ref
# idk what that means but it wasn't an error soooooooo
design_bottom <- modelMatrix(targets_bottom, ref="Ref")
# > Found unique target names:
  # > Cancer PSC Ref

## estimate array quality weights: ----
# arrayWeights "estimates relative quality weights"
aw_top <- arrayWeights(MAtopx, design_top)
aw_bottom <- arrayWeights(MAbottomx, design_bottom)

## account for duplicate spots ----
# duplicateCorrelation "estimates intra-block correlation"
# ndups is the number of duplicates on the slides
dupcor_top <- duplicateCorrelation(MAtopx, design_top, ndups = 2)
dupcor_bottom <- duplicateCorrelation(MAbottomx, design_bottom, ndups = 2)
# shrug emoji

## fit a linear model ----
# "fits a linear model for each gene given a series of arrays"
fit_top <- lmFit(MAtopx, design = design_top, ndups = 2, correlation = dupcor_top$consensus.correlation)
# > Warning message: Partial NA coefficients for 54 probe(s)
# it's a warning not an error sooooooooo
fit_bottom <- lmFit(MAbottomx, design = design_bottom, ndups = 2, correlation = dupcor_bottom$consensus.correlation)
# > Partial NA coefficients for 5 probe(s) 

## make a contrast matrix ----
# CHANGE : the samples i want to compare are (in this case) PSC & Cancer, labelled as such in the Cy5 column in the targets files
contrast.matrix_top <- makeContrasts(CRC_sc-PDAC_sc, levels = design_top)
contrast.matrix_bottom <- makeContrasts(CRC_sc-PDAC_sc, levels = design_bottom)

contrast.matrix_top <- makeContrasts(HCC_sr-PDAC_sr, levels = design_top)
contrast.matrix_bottom <- makeContrasts(HCC_sr-PDAC_sr, levels = design_bottom)


## now combining the fit and the contrasts ----
fit_top1 <- contrasts.fit(fit_top, contrast.matrix_top)
fit_bottom1 <- contrasts.fit(fit_bottom, contrast.matrix_bottom)
# well there are no errors. onwards

## application of ebayes ----
# computes moderated t-stats, f-stats, and log-odds of differential exp by empirical bayes moderation of the standard errors towards a global value
fit_top2 <- eBayes(fit_top1)
# topTable extracts a table of the top ranked genes from a linear model fit
dif_top <- topTable(fit_top2, coef = 1, number = nrow(fit_top2), adjust.method = "BH")
# coef is the column of interest, n = max # of genes to list, adjust.method is BH (benjamini & hochberg)
#repeating for the bottom
fit_bottom2 <- eBayes(fit_bottom1)
# topTable extracts a table of the top ranked genes from a linear model fit
dif_bottom <- topTable(fit_bottom2, coef = 1, n = nrow(fit_bottom2), adjust.method = "BH")

#in viewing all the things, idk what the numbers are sorted by, bc it's definitely not p-value or adjusted p-value

## volcano plot time! ----
volcanoplot(fit_top2, coef = 1, highlight = 10, names = fit_top2$genes$Name, main = "Volcano Plot Top") # => volcano_plot_top
# could do with playing with the plot layout but hey results are there
volcanoplot(fit_bottom2, coef = 1, highlight = 10, names = fit_bottom2$genes$Name, main = "Volcano Plot Bottom") # => volcano_plot_bottom
# again, there exists results


## exporting the data, hopefully it means something ----
write.csv(dif_top, file = "CRCsc-PDACsc_top.csv")
write.csv(dif_bottom, file = "CRCsc-PDACsc_bottom.csv")


write.csv(dif_top, file = "HCCsr-PDACsr_top.csv")
write.csv(dif_bottom, file = "HCCsr-PDACsr_bottom.csv")
## so now the actual results??? ----
# decideTests -> identifies which genes are statistically differentially exp for each contrast
results_top <- decideTests(fit_top2, method = "global")
results_bottom <- decideTests(fit_bottom2, method = "global")
# not sure how we use these results bc it just looks like a table w -1s, 0s, and 1s
vennDiagram(results_top, include = "up", cex.main = 2, main = "Top - upreg", line = -3) # => top_upreg_venn
vennDiagram(results_top, include = "down", main = "TOP - downreg", cex.main = 2, line = -3) # => top_downreg_venn
vennDiagram(results_bottom, include = "up", cex.main = 2, main = "Bottom - upreg", line = -3) # => bottom_upreg_venn
vennDiagram(results_bottom, include = "down", main = "Bottom - downreg", cex.main = 2, line = -3) # => bottom_downreg_venn
# loads of warning messages but the diagrams still show up soooo

# that's that i guess? gonna talk to fawaz abt the results lol