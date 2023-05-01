# Modified standardizedPosterior() function
# Modification is that it allows the draws to be replaced

standardizedPosterior.replace.draws <- function (object, replacement.draws, ...)
{
  dots <- list(...)
  allowargs <- c("type", "cov.std", "remove.eq",
                 "remove.ineq", "remove.def")
  if (!all(names(dots) %in% allowargs)) {
    stop(paste0("blavaan ERROR: arguments must be in ",
                paste(allowargs, collapse = " ")))
  }

  draws.from.bsem.object <- make_mcmc(object@external$mcmcout)
  draws.from.bsem.object <- do.call("rbind", draws.from.bsem.object)
  n.draws.from.bsem <- nrow(draws.from.bsem.object)

  replaced.draws <- rbind(
    draws.from.bsem.object,
    matrix(NA, nrow=nrow(replacement.draws)-n.draws.from.bsem, ncol=ncol(draws.from.bsem.object))
  )

  for(which.col in 1:ncol(replaced.draws)){
    column.name <- colnames(replaced.draws)[which.col]

    # check if the column name in the bsem object is also one for which there are draws from the combined model from jags
    # and if so, replace it
    # Otherwise, extend down the values of that column based on what's in the first row
    if(any(column.name == colnames(replacement.draws))){
      replaced.draws[ , column.name] <- replacement.draws[ , column.name]
    }
    if(!any(column.name == colnames(replacement.draws))){
      replaced.draws[ , column.name] <- replaced.draws[1, column.name]
    }
    #which(column.name == colnames(measurement.model.draws.as.data.frame))
  }


  #draws <- make_mcmc(object@external$mcmcout)
  #draws <- do.call("rbind", draws)
  draws <- replaced.draws

  tf <- object
  tmp <- fill_params(draws[1, ], object@Model, object@ParTable)
  tf@Model <- tmp
  tf@ParTable$est[tf@ParTable$free > 0] <- lav_model_get_parameters(tmp)
  tmp2 <- do.call("standardizedSolution", c(list(object = tf),
                                            dots))
  fullres <- matrix(NA, nrow(draws), nrow(tmp2))
  colnames(fullres) <- with(tmp2, paste0(lhs, op, rhs))
  if ("group" %in% colnames(tmp2))
    colnames(fullres) <- paste0(colnames(fullres), ".g",
                                tmp2$group)
  fullres[1, ] <- tmp2[, "est.std"]
  for (i in 2:nrow(draws)) {
    tmp <- fill_params(draws[i, ], object@Model, object@ParTable)
    tf@Model <- tmp
    tf@ParTable$est[tf@ParTable$free > 0] <- lav_model_get_parameters(tmp)
    fullres[i, ] <- do.call("standardizedSolution",
                            c(list(object = tf), dots))[, "est.std"]
  }
  fullres
}



