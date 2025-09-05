delta_metric <- function(cor_matrix){

  # cor_matrix contains correlations for different cell types (rows) at different metric settings (columns)
  # For example:
  #                  0   0.01   0.05    
  # Endothelium  0.541  0.581  0.167  
  # Epithelial   0.756  0.725  0.684 
  # Muscle       0.881  0.869  0.780
  # 
  # Here columns are different lambda settings and rows are different cell types.
  
  nc = ncol(cor_matrix)
  nr = nrow(cor_matrix)
  
  delta1 = matrix(NA, nrow= nr, ncol=nc)
  delta2 = matrix(NA, nrow= nr, ncol=nc)
  celltype_means1 = as.numeric(rowMeans(cor_matrix)) # regular mean
  
  for (i in 1:nr){
    celltype_means2 = mean(sort(as.numeric(cor_matrix[i,]))[2:(nc-1)]) # winsorized mean
    for(j in 1:nc){
      delta1[i,j] = cor_matrix[i,j] - celltype_means1[i]
      delta2[i,j] = cor_matrix[i,j] - celltype_means2
    }
  }
  
  delta_lamda1 = colMeans(delta1) 
  delta_lamda2 = colMeans(delta2)
  delta_lamda3 = rep(NA, length = length(delta_lamda1))
  delta_lamda4 = rep(NA, length = length(delta_lamda1))
  for(i in 1:ncol(delta1)){
    temp = delta1[,i]
    temp[temp>0.25]=0.25
    temp[temp<-0.25]=-0.25
    delta_lamda3[i]=mean(temp)
    temp = delta2[,i]
    temp[temp>0.25]=0.25
    temp[temp<-0.25]=-0.25
    delta_lamda4[i]=mean(temp)
  }
  
  results = data.frame(metrics = colnames(cor_matrix),
                       delta = round(delta_lamda1,digits= 4),
                       delta_robust = round(delta_lamda2, digits = 4),
                       delta_winsorized = round(delta_lamda3,digits= 4),
                       delta_robust_winsorized = round(delta_lamda4,digits= 4))
  
  return(results)
}
