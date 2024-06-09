library(spacexr)
library(Matrix)
library(DescTools)
counts <- read.csv("CerebellumPuck_sc_count.tsv",
                   sep = "\t", stringsAsFactors=FALSE)# load in counts matrix
rownames(counts) <- counts[,1]
counts[,1] <- NULL
counts <- t(counts)
meta_data <- read.csv("CerebellumPuck_sc_label.tsv",
                   sep = "\t")
cell_types <- meta_data$bio_celltype
names(cell_types) <- meta_data$cell # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
### Create the Reference object
reference <- Reference(counts, cell_types)
#saveRDS(reference,"reference_mouse.rds")
## Load ST data
counts <- read.csv("Cerebellum_counts.csv") # load in counts matrix
coords <- read.csv("Cerebellum_coords.csv")
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords$barcode; coords$barcode <- NULL # Move barcodes to rownames
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)
myRCTD <- create.RCTD(puck, reference, max_cores = 4, 
                      UMI_min = 1,UMI_max=max(nUMI))
# full mode
full <- run.RCTD(myRCTD, doublet_mode = 'full') 
results <- full@results
normalized_weights <- normalize_weights(results$weights)
H_mod=t(normalized_weights) # estimated proportions for cell types
