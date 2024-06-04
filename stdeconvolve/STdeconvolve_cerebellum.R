# record memory usage
gc()
Rprof("Rprof.out", memory.profiling=TRUE)
#record time
start.time <- Sys.time()
#####################################
library(STdeconvolve)
library(corrplot)
library(dplyr)
library(ggplot2)
library(cowplot)
library(DescTools)
library(retrofit)
setwd("/storage/group/qul12/default/retrofit/Mouse/")
path <- "Cerebellum_counts.csv"
ppos <- "Cerebellum_coords.csv"
pos <- read.csv(file = ppos, row.names = 1)
real <- read.csv(file = path, row.names = 1)
#input: pixel by gene
spot <- cleanCounts(counts = as.matrix(real),plot = TRUE,verbose = TRUE)
# 61 genes
gene <- restrictCorpus(spot,plot = TRUE,verbose = TRUE)
counts <- t(as.matrix(gene))
mod <- fitLDA(counts[rowSums(counts)>0,], 
              #as.matrix(real_all),
              Ks = 20, 
              #Ks =seq(2,20,by=1),
              plot=TRUE, verbose=TRUE)
mod30 <- optimalModel(models = mod, opt = 20)
res30 <- getBetaTheta(mod30, betaScale = 1000)
######################################
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
Rprof(NULL)
peak <- max(summaryRprof("Rprof.out", memory="both")$by.total$mem.total)
print(paste0(peak, "Mb"))