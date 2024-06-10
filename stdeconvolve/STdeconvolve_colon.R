rm(list=ls())
library(STdeconvolve)
library(corrplot)
library(dplyr)
library(ggplot2)
library(cowplot)
library(DescTools)
library(retrofit)

pos <- read.csv("A3_coords.csv", row.names = 1)
real <- read.csv("A3_counts.csv", row.names = 1)
#input: pixel by gene
spot <- cleanCounts(counts = as.matrix(real),plot = TRUE,verbose = TRUE)
# 61 genes
gene <- restrictCorpus(spot,plot = TRUE,verbose = TRUE)
counts <- t(as.matrix(gene))
mod <- fitLDA(counts[rowSums(counts)>0,], 
              #as.matrix(real_all),
              Ks = 16, 
              #Ks =seq(2,20,by=1),
              plot=TRUE, verbose=TRUE)
mod <- optimalModel(models = mod, opt = 16)
res <- getBetaTheta(mod, betaScale = 1000)
rmProp <- res$theta # estimated proportions for cell types
rmGexp <- res$beta # estimated gene expression
