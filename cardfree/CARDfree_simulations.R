library(CARD)
library(retrofit)
library(DescTools)
library(corrplot)
##load ST and location files
loc <- "BeadLocationsForR.csv"
location <- read.csv(file = loc,row.names = 1)
colnames(location) <- c("x", "y")

input_x = read.csv("N=10,M=3_loc_X.csv")
rownames(input_x) = input_x[,1]
input_x = as.matrix(input_x[,-1])
# Create marker gene list: od gene + sc expression
gt_W0 = read.csv("Cerebellum_W_K=10.csv")
rownames(gt_W0) = gt_W0[,1]
gt_W0 = as.matrix(gt_W0[,-1])
de_gene = read.table("187od_gene_N=10,M=3.txt")
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

###Annotation
gt_w0_sub <- gt_W0[unlist(de_gene),]
gt_w0_sub <- gt_w0_sub[order(row.names(gt_w0_sub)),]
## match the latent components to known cell types using a cell-type-specific gene expression reference
res <- retrofit::annotateWithCorrelations(gt_w0_sub, K=10, 
            decomp_w=CARDfree_obj@estimated_refMatrix, #decomp_w Matrix(GeneExpressions, Components): Decomposed w matrix
            decomp_h=t(CARDfree_obj@Proportion_CARD) #decomp_h Matrix(Components, Spots): Decomposed h matrix 
)
W_prop <- res$w_prop # estimated gene expression profiles
H_prop <- res$h_prop # estimated cell type proprotions
