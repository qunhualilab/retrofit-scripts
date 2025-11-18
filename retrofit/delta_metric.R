delta_metric <- function(cor_matrix){
  
  #' RETROFIT delta metric 
  #' 
  #' @description Receiving the input with correlation matrix for different cell types (rows) at different parameter settings (columns), 
  #'  the function returns the delta metric.
  #'
  #' @param cor_matrix matrix or array with dimension (CellTypes, Parameters). 
  #' For example: cor_matrix =
  #'                  0   0.01   0.05    
  #' Endothelium  0.541  0.581  0.167  
  #' Epithelial   0.756  0.725  0.684 
  #' Muscle       0.881  0.869  0.780
  #' 
  #' Here, there are 3 cell types and 3 different parameter settings. 
  #' The values of the matrix are the correlation between the estimated gene expression profile via RETROFIT,
  #' and the corresponding cell-type specific marker gene expressions at different parameter settings.
  #'
  #' @return A dataframe that contains
  #' \itemize{
  #' \item parameters: different parameter settings provided in the input
  #' \item delta:      delta metric for each parameter setting
  #' }
  #' 
  #'@seealso Algorithm 1 in paper
  #'@export
  
  nc = ncol(cor_matrix)
  nr = nrow(cor_matrix)
  
  delta = matrix(NA, nrow= nr, ncol=nc)
  celltype_means = as.numeric(rowMeans(cor_matrix)) 
  
  for (i in 1:nr){
    for(j in 1:nc){
      delta[i,j] = cor_matrix[i,j] - celltype_means[i]
    }
  }
  
  delta_param = colMeans(delta) 
  
  results = data.frame(parameters = colnames(cor_matrix),
                       delta = round(delta_param,digits= 4))
  
  return(results)
}
