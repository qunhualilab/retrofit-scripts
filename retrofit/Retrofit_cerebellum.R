## RETROFIT for Mouse Brain Colon data ##
####################################
rm(list=ls())
library(retrofit)

## Load ST data 
X = read.csv("Cerebellum_X.csv", row.names = 1)

## Load references
sc_ref = read.csv("Cerebellum_all_W.csv", row.names = 1) # single-cell reference
df=read.csv("CerebellumPuck_markers.csv", row.names = 1) # marker reference

## fixing gene names in SC reference
for(i in 1:nrow(sc_ref)){
  if(rownames(sc_ref)[i]=="Hbb-bs"){
    rownames(sc_ref)[i] ="Hbb.bs"
  }
  if(rownames(sc_ref)[i]=="mt-Rnr2"){
    rownames(sc_ref)[i] ="mt.Rnr2"
  }
}
sc_ref = sc_ref[rownames(sc_ref) %in% rownames(X),] # keeping only common genes in reference
colnames(sc_ref)=c("Astrocytes","Bergmann Glia", "Choroid" ,
                   "Endothelial", "Granule" ,"Interneurons","Microglia" ,
                   "Muraland Tip"  ,"Oligo","Purkinje",
                   "PV Interneurons") 

## creating marker reference
marker_ref <- list(Purkinje = df$Gene[df$Celltype=="Purkinje"],
                   Endothelial = df$Gene[df$Celltype=="Endothelial"],
                   Oligo = df$Gene[df$Celltype=="Oligo"],
                   Bergmann_Glia = df$Gene[df$Celltype=="Bergmann Glia"],
                   Choroid = df$Gene[df$Celltype=="Choroid"],
                   Granular = df$Gene[df$Celltype=="Granule"],
                   Interneurons = df$Gene[df$Celltype=="Interneurons"],
                   Muraland_Tip = df$Gene[df$Celltype=="Muraland Tip"],
                   Astrocytes = df$Gene[df$Celltype=="Astrocytes"],
                   Microglia = df$Gene[df$Celltype=="Microglia"],
                   PV_Interneurons = df$Gene[df$Celltype=="PV Interneurons"])




## Setting parameters
L=20
iter = 4000
lamda = 0.05

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

