# Residual Correlations for all pairings of variables in a matrix
# Arguments: 
#   data.matrix: data matrix
#   mod.imp.cov.matrix: model-implied covariance matrix
#   output: character taking on possibilities:
#     "cor matrix" 
#     "long"        is long format, with 'var1' and 'var2' columns indicating the pairings of observables, and the 'cor' column indicating the residual correlation
#     "wide"        is wide format, with 1 row, and each column is pairing of variables (e.g., x1_x2)
#     default is "wide"

#data.matrix = observables.to.analyze
#mod.imp.cov.matrix = mod.imp.cov.matrix.observables
#output="wide"

resid.cor.matrix.function <- function(data.matrix, mod.imp.cov.matrix, output="wide"){

  # compute model-implied standard deivations
  # and store them on the diagonal of a matrix
  mod.imp.sds <- sqrt(diag(mod.imp.cov.matrix))
  D <- diag(ncol(mod.imp.cov.matrix))
  diag(D) <- mod.imp.sds
  
  # compute model-implied correlation matrix
  mod.imp.cor.matrix <- solve(D) %*% mod.imp.cov.matrix %*% solve(D)
  
  # Calculate the residual correlation matrix
  resid.cor.matrix = cor(data.matrix) - mod.imp.cor.matrix
  
  if(output=="cor matrix"){
    result=resid.cor.matrix
  }
  
  if(output=="long"){
    result=resid.cor.matrix %>% rstatix::cor_gather()
  }
  
  if(output=="wide"){
    result=resid.cor.matrix %>% rstatix::cor_gather()
    pairing.name <- paste0(result$var1, "_", result$var2)
    result <- cbind(result, pairing.name)
    result <- result %>% tidyr::pivot_wider(names_from = pairing.name, values_from = cor)
    result <- select(result, -var1)
    result <- select(result, -var2)
    result <- colMeans(result, na.rm=TRUE)
    result
  }
    
  result
}