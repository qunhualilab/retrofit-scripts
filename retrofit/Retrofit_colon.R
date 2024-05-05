## RETROFIT for Human Colon data ##
####################################
rm(list=ls())
library(retrofit)

## Load ST data 
X = read.csv("A3_X.csv", row.names = 1)

## Load references
sc_ref = read.csv("Intestine_W_12pcw.csv", row.names = 1) # single cell reference
sc_ref = sc_ref[rownames(sc_ref) %in% rownames(X),] # keeping only common genes in reference
celltypes = c('Endothelium','Epithelial','Fibroblasts','Immune','Muscle','Myo.Meso',
              'Neural','Pericytes')
sc_ref = sc_ref[,celltypes] # fixing cell type order

df = read.csv("Colon_markers.csv")[,-1] # marker reference
marker_ref <- list(Endothelium = df$Gene[df$Celltype=="Endothelium"],
                   Epithelial = df$Gene[df$Celltype=="Epithelial"],
                   Fibroblasts = df$Gene[df$Celltype=="Fibroblasts"],
                   Immune = df$Gene[df$Celltype=="Immune"],
                   Muscle = df$Gene[df$Celltype=="Muscle"],
                   Myo.Meso = df$Gene[df$Celltype=="Myo.Meso"],
                   Neural = df$Gene[df$Celltype=="Neural"],
                   Pericytes = df$Gene[df$Celltype=="Pericytes"])

## Setting hyper parameters
L=16
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
sc_ref = sc_ref[rownames(W),]
result1    = annotateWithCorrelations(sc_ref, K, W, H)
H_sc_prop = result1$h_prop # estimated proportions for cell types
W_sc_prop = result1$w_prop # estimated gene specific information

## Mapping with markers
result2      <- annotateWithMarkers(marker_ref=marker_ref, K, W, H)
H_marker_prop = result2$h_prop # estimated proportions for cell types
W_marker_prop = result2$w_prop # estimated gene specific information