rm(list=ls())
library(STdeconvolve)
library(corrplot)
library(dplyr)
library(ggplot2)
library(cowplot)
library(DescTools)
library(retrofit)

pos <- read.csv("BeadLocationsForR.csv", row.names = 1)
real <- read.csv("N=10,M=3_loc_X.csv", row.names = 1)
#input: pixel by gene
spot <- cleanCounts(counts = as.matrix(real),plot = TRUE,verbose = TRUE)
# 187 genes
gene <- restrictCorpus(spot,plot = TRUE,verbose = TRUE)
counts <- t(as.matrix(gene))
mod <- fitLDA(counts[rowSums(counts)>0,], 
              #as.matrix(real_all),
              Ks = 20, 
              #Ks =seq(2,20,by=1),
              plot=TRUE, verbose=TRUE)
modres <- optimalModel(models = mod, opt = 20)
res10 <- getBetaTheta(modres, betaScale = 1000)
rmProp <- res10$theta # estimated proportions for cell types
rmGexp <- res10$beta # estimated gene expression
# annotation
ptruexp <- "Cerebellum_W_K=10.csv"
truexp <- read.csv(file = ptruexp, row.names = 1)
genelist <- colnames(rmGexp)
te <- truexp[genelist,]
H = t(rmProp)
W = t(rmGexp)

## match the latent components to known cell types using a cell-type-specific gene expression reference
res1 <- retrofit::annotateWithCorrelations(te, K=10, decomp_w=W, decomp_h=H)
h_1sc <- res1$h
H_1sc_prop = res1$h_prop # estimated cell type proportion
W_1sc_prop = res1$w_prop # estimated gene expression
