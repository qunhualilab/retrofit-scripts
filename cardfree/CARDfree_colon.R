library(CARD)
library(retrofit)
library(DescTools)
library(corrplot)
##load ST and location files
loc <- "A3_coords.csv"
location <- read.csv(file = loc,row.names = 1)
colnames(location) <- c("x", "y")

input_x = read.csv("A3_counts.csv")
rownames(input_x) = input_x[,1]
input_x = as.matrix(input_x[,-1])
# Create marker gene list: od gene + sc expression
marker <- "Colon_Markers.csv"
df=read.csv(marker, row.names = 1)
df$Celltype[which(df$Celltype=="Endothelium")]<- "Endothelial"
df$Celltype[which(df$Celltype=="Fibroblasts")]<- "Fibroblast"
df$Celltype[which(df$Celltype=="Myo.Meso")]<- "MyoFB/Meso"
df$Celltype[which(df$Celltype=="Pericytes")]<- "Pericyte"

gt_W0 = read.csv("Intestine_W_12pcw.csv")
rownames(gt_W0) = gt_W0[,1]
gt_W0 = as.matrix(gt_W0[,-1])# [,c(1:5)]

de_gene = read.table("od_gene_ColonA3.txt") 
de_gene_m = c(de_gene[,1],df$Gene)
sub_W = gt_W0[rownames(gt_W0) %in% de_gene_m,]
cell_type = apply(sub_W, 1, which.max)

mk_gene = list()
for(i in 1:ncol(gt_W0)){
  print(colnames(gt_W0)[i])
  mk_gene = c(mk_gene, list(names(which(cell_type == i))))
}
names(mk_gene) = colnames(gt_W0)

CARDfree_obj = createCARDfreeObject(markerList = mk_gene, 
                                    spatial_count = input_x,  
                                    spatial_location = location,
                                    minCountGene = 0,
                                    minCountSpot = 0)

## deconvolution using CARDfree
CARDfree_obj = CARD_refFree(CARDfree_obj)
print(CARDfree_obj@Proportion_CARD[1:2,])
