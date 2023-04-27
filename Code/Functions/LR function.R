# Likelihood ratio discrepancy measure for a covariance matrix
LR.function <- function(data.cov.matrix, mod.imp.cov.matrix, n){

  F.ML <- log(det(mod.imp.cov.matrix)) + sum(diag((data.cov.matrix %*% solve(mod.imp.cov.matrix)))) - log(det(data.cov.matrix)) - ncol(data.cov.matrix)
  LR <- (n-1)*F.ML
	LR

}