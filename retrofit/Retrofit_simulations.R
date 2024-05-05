## RETROFIT for simulated data ##
####################################
rm(list=ls())
library(retrofit)

## Load ST data 
X = read.csv("N=10,M=3_loc_X.csv", row.names = 1)

## Load references
sc_ref = read.csv("Cerebellum_W_K=10.csv", row.names = 1)

## Setting parameters
L=20
G = nrow(X) # number of genes
S = ncol(X) # number of spots
iter = 4000
lamda = 0.01

## Deconvolution
set.seed(1)
result = retrofit::decompose(X, L=L, iterations=iter, lambda = lamda, verbose=FALSE)

H   = result$h
W   = result$w
Theta = result$th

## Mapping using single cell
K = ncol(sc_ref) # number of cell types
result    = annotateWithCorrelations(sc_ref, K, W, H)
H_sc_prop = result$h_prop # estimated proportions for cell types
W_sc_prop = result$w_prop # estimated gene specific information

