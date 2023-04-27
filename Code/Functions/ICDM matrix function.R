# Interpretational confouding discrepancy measures for all pairings of variables in a matrix
# Arguments: 
#   indicators.data.matrix: data matrix of the indicators
#   indicators.loadings: vector of loadings of indicators on (single) latent variable
#   non.indicators.data.matrix: data matrix of the non-indicators, typically covariates or outcomes

#indicators.data.matrix = indicators
#if(model.has.covariates) non.indicators.data.matrix = covariates
#if(model.has.outcomes) non.indicators.data.matrix = outcomes
#indicators.loadings = loadings.this.iter

ICDM.function <- function(indicators.data.matrix, indicators.loadings, non.indicators.data.matrix){
  
  # set up objects to collect values and names
  result=NULL
  names=NULL
  
  # find number of indicators, external variables for each pair of indicators to be paired with
  J.local=ncol(indicators.data.matrix)
  num.external.vars <- ncol(non.indicators.data.matrix)
  
  # loop over triplets of one indicator, other indicator, and external variable
  for(j in 1:(J.local-1)){
    for(jj in (j+1):J.local){
      for(which.external.var in 1:num.external.vars){
        
        # compute ICDM
        temp = (cov(non.indicators.data.matrix[,which.external.var], indicators.data.matrix[,j])/cov(non.indicators.data.matrix[,which.external.var], indicators.data.matrix[,jj])) - (indicators.loadings[j]/indicators.loadings[jj])
        temp = as.numeric(temp)
        result=cbind(result, temp)

        # Define name
        names <- c(names,paste0(
          "ICDM_", 
          colnames(indicators.data.matrix)[j],
          "_",
          colnames(indicators.data.matrix)[jj],
          "_",
          colnames(non.indicators.data.matrix)[which.external.var]
        )
        )
        
      } # closes loop over which.external.var
    } # closes loop over jj
  } # closes loop over j
  
  colnames(result) <- names
  result
}