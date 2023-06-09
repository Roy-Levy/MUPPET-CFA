# Governing R code to run analyses for
# Measurement and Uncertainty Preserving Parametric (MUPPET) CFA on
# The IIS Example
# Using the functions


# 4/25/2023
# Testing how to hopefully use blavaan to automate use of JAGS in fitting structural portion
# Modified version of 'IIS Example Using Functions (mine).R' from 1/26/2022


rm(list=ls())



# Define the folders ------

#   Drive letter
#   main analysis folder
#   R Code folder and functions subfolder
#   Conditions folder

main.folder <- getwd()

# library(rprojroot)
# file.name.in.root.directory <- "MUPPET CFA IIS (Testing).Rproj"
# main.folder <- find_root(has_file(file.name.in.root.directory))

code.folder <- paste(main.folder, "/Code/", sep="")
functions.folder <- paste(code.folder, "Functions/", sep="")

data.folder <- paste(main.folder, "/Data/", sep="")


# Read in and prepare the data only need to do this once ------
if(1==0){
# * Read in the data
file.name <- "IIS Final Modality Ground 10-16-08.sav"
raw.data <- haven::read_spss(file=paste0(data.folder, file.name))


# * Prepare the data
library(dplyr)

working.data <- raw.data
working.data <- working.data %>% select(BIRTHDATE, V1_mean, V2_mean, V3_mean, V4_mean, V5_mean)

# remove extra text from birthdate
working.data$BIRTHDATE <- stringr::str_remove(working.data$BIRTHDATE, " 00:00:00")

# Convert birthdate to date forma
working.data$BIRTHDATE <- as.Date(working.data$BIRTHDATE, "%m/%d/%Y")

# Define the date at which data were collected
dataset.data <- as.Date("2008-10-16")

# Add age to the dataset
working.data <- mutate(
  working.data,
  age_in_years=as.numeric(
    difftime(time1 = dataset.data, working.data$BIRTHDATE, units = "weeks"))*0.019165
  )

# Rename variables
colnames(working.data) <- c("birthdate", "PI", "FI", "AD", "FC", "IGC", "age_in_years")


# Declare the variables that are indicators, covariates, and outcomes ------
indicators.names <- c("PI", "FI", "AD", "FC")
covariates.names <- "age_in_years"
outcomes.names <- "IGC"


# Define the data to analyze, including order of variables ------

# Define the observables
all.observables <- select(working.data, contains(c(covariates.names, indicators.names, outcomes.names)))

# Reorder the observables
all.observables <- select(all.observables, c(covariates.names, indicators.names, outcomes.names))

# Remove missing cases
all.observables.complete.data <- na.omit(all.observables)

# Write out the data ---------
file.name <- "IIS.dat"
write.table(all.observables.complete.data, paste0(data.folder, file.name), row.names = FALSE)
} # closes switch for prepping data, only need to do once



# Load the packages needed -----
# packages
library(dplyr)
library(blavaan)
future::plan("multisession")
options(mc.cores = parallel::detectCores())

library(MCMCvis)
library(coda)

library(parallel)
library(foreach)
library(doParallel)

library(MASS)
library(lavaan)

library(rstatix)

library(stringr)


# Load the functions -----
file.name <- "LR function.R"
source(paste0(functions.folder, file.name))

file.name <- "Residual correlation matrix function.R"
source(paste0(functions.folder, file.name))

file.name <- "ICDM matrix function.R"
source(paste0(functions.folder, file.name))

file.name <- "standardizedPosterior.replace.draws function.R"
source(paste0(functions.folder, file.name))

file.name <- "blav_utils.R"
source(paste0(functions.folder, file.name))

file.name <- "MUPPET CFA function (mine).R"
source(paste0(functions.folder, file.name))


# Read in the data ------
file.name <- "IIS.dat"
raw.data <- as.data.frame(read.table(paste0(data.folder, file.name), header = TRUE))


data=raw.data
indicators.names=c("PI", "FI", "AD", "FC")
covariates.names=NULL
outcomes.names="IGC"
measurement.model.priors=dpriors(
  #nu="normal(0,3)",        # measured variable intercepts in N(mean, sd)
  lambda="normal(0,10)",    # measured variable loadings in N(mean, sd)
  #theta="gamma(.01,.01)[prec]"     # measured variable error precision
  theta="gamma(1,1)[prec]"     # measured variable error precision
  #psi="gamma(5,10)[prec]"        # latent variable precision
)
n.chains = 2
n.warmup = 500
n.burnin = 0
n.iters.per.chain.after.warmup.and.burnin = 10
obtain.standardized.cfa=TRUE
obtain.standardized.combined.model = TRUE
beta.o.prior.mean = 0
beta.o.prior.var = 10000
psi.y.prior.alpha = 1
psi.y.prior.beta = 1
beta.c.prior.mean = 0
beta.c.prior.var = 10000
n.iters.per.chain.total.structural = 51
save.summary.stats.from.MUPPET=TRUE
save.draws.from.MUPPET=TRUE
model.check=FALSE
save.post.pred.data=FALSE


# First fit the CFA in stage 1, here not using a function ------
# So just declare values of arguments
data=raw.data
indicators.names=c("PI", "FI", "AD", "FC")
covariates.names=NULL
outcomes.names="IGC"
measurement.model.priors=dpriors(
  #nu="normal(0,3)",        # measured variable intercepts in N(mean, sd)
  lambda="normal(0,10)",    # measured variable loadings in N(mean, sd)
  #theta="gamma(.01,.01)[prec]"     # measured variable error precision
  theta="gamma(1,1)[prec]"     # measured variable error precision
  #psi="gamma(5,10)[prec]"        # latent variable precision
)
n.chains = 2
n.warmup = 500
n.burnin = 0
n.iters.per.chain.after.warmup.and.burnin = 5
obtain.standardized.cfa=TRUE
beta.o.prior.mean = 0
beta.o.prior.var = 10000
psi.y.prior.alpha = 1
psi.y.prior.beta = 1
beta.c.prior.mean = 0
beta.c.prior.var = 10000
n.iters.per.chain.total.structural = 51
save.summary.stats.from.MUPPET=TRUE
save.draws.from.MUPPET=TRUE
model.check=FALSE
save.post.pred.data=FALSE


# Below is ripped from the function
# * Center the data ------
data.centered <- as.data.frame(scale(data, scale=FALSE))

# * Declare if covariates or outcomes are present ------
model.has.covariates <- length(covariates.names)!=0
model.has.outcomes <- length(outcomes.names)!=0

# * Datasets for indicators, covariates, and outcomes  ------
indicators <- as.data.frame(select(data.centered, all_of(indicators.names)))
J = length(indicators.names)
colnames(indicators) <- paste0("x", seq(1:J))
if(model.has.covariates) {
  covariates <- as.data.frame(dplyr::select(data.centered, all_of(covariates.names)))
  H = length(covariates.names)
  colnames(covariates) <- paste0("w", seq(1:H))
}
if(model.has.outcomes) {
  outcomes <- as.data.frame(dplyr::select(data.centered, all_of(outcomes.names)))
  K = length(outcomes.names)
  colnames(outcomes) <- paste0("y", seq(1:K))
}

# * Declare number  latent variables ------
M=1

# STAGE 1: MEASUREMENT MODEL ---
# Load packages for fitting the measurement model ----
#library(blavaan)
#future::plan("multisession")
#options(mc.cores = parallel::detectCores())

# Generate lavaan syntax ------

# * Measurement model factor structure -----
cfa.model.syntax.lavaan <- 'f1 =~ NA*'
cfa.model.syntax.lavaan <- paste0(cfa.model.syntax.lavaan, colnames(indicators)[1])
if(J > 1){
  for(j in 2:J){
    cfa.model.syntax.lavaan <- paste0(cfa.model.syntax.lavaan, " + NA*", colnames(indicators)[j])
  }
} # closes if J>1

# * Fix intercepts to 0 in measurement model  -----
for(j in 1:J){
  cfa.model.syntax.lavaan <- paste(cfa.model.syntax.lavaan, paste0(colnames(indicators)[j], " ~ 0*1"), sep="\n")
}

# cat(cfa.model.syntax.lavaan)

# Fit the measurement model ------

n.iters.total.per.chain = n.warmup+n.burnin+n.iters.per.chain.after.warmup.and.burnin

# Print out beginning of fitting this portion of the model
print(paste0("Fitting the measurement model for ",
             n.chains, " chains and ",
             n.iters.total.per.chain, " total iterations per chain")
)


fitted.model.bcfa <- bcfa(
  model = cfa.model.syntax.lavaan,
  dp=measurement.model.priors,
  n.chains=n.chains,
  #burnin = n.burnin,
  burnin = n.warmup,   # Think this is warmup in stan
  #adapt=n.warmup,
  sample=n.iters.per.chain.after.warmup.and.burnin,
  std.lv = TRUE,			# identify the model by fixing factor variance to 1
  #std.ov = TRUE,      # to get standardized solution
  int.ov.free = FALSE,
  #estimator = "ML",
  #likelihood = "wishart", # normal based  on dividing by N, wishart is N-1
  #se = "standard"
  mcmcfile=TRUE,
  save.lvs=FALSE,
  test="none",   # turn off computing fit functions
  bcontrol=list(cores=n.chains),
  data = indicators
)


# Obtain the standardized solution to the measurement model -----
# if(obtain.standardized.cfa) standardized.posterior.cfa <- standardizedPosterior(fitted.model.bcfa)


# STAGE 2: STRUCTURAL MODEL ----

# Load the packages for preparing draws from previous stage for next stage -----
#library(MCMCvis)
#library(coda)
#library(dplyr)

# Extract the stanfit object -----
fitted.model.stanfit <- blavInspect(fitted.model.bcfa, what="mcobj")

# convert draws to mcmc.list
draws.as.mcmc.list <- MCMCchains(fitted.model.stanfit,
                                 mcmc.list = TRUE)

# Define the draws to analyze
draws.to.analyze <- draws.as.mcmc.list

# Extract the parameter table from blavaan for the CFA model using stan
cfa.partable.stan <- as.data.frame(fitted.model.bcfa@ParTable)

# Define the internal names of the estimated parameters in stan
parameter.name.pxnames <- (cfa.partable.stan %>% filter(!is.na(pxnames)))$pxnames


# Convert the draws from the measurement model to a data frame ------
# measurement.model.draws.as.data.frame <- dplyr::select(
#   as.data.frame(as.matrix(draws.to.analyze)),
#   contains(c("ly","Theta"))
# )

measurement.model.draws.as.data.frame <- dplyr::select(
  as.data.frame(as.matrix(draws.to.analyze)),
  contains(parameter.name.pxnames)
)




# Get what we think blavaan will call the parameters in JAGS, from the stan results -----
parameter.name.jags <- c()
temp <- cfa.partable.stan %>% filter(!is.na(stanpnum))

for(which.row in 1:nrow(temp)){
  parameter.name.jags <- c(parameter.name.jags,
                           paste0(
                             temp[which.row,]$mat,
                             "[",
                             temp[which.row,]$row,
                             ",",
                             temp[which.row,]$col,
                             ",",
                             temp[which.row,]$group,
                             "]"
                           )
  )
}

cfa.partable.stan.estimated.parameters <- cbind(temp, parameter.name.jags)


### GO TO SCRATCH FILE -----



# Add on the values for the factor variance of 1
phi <- rep(1, nrow(measurement.model.draws.as.data.frame))
measurement.model.draws.as.data.frame <- cbind(measurement.model.draws.as.data.frame, phi)



# Try to fix parameters in JAGS ------------------
# testing fitting the model in JAGS, just to see JAGS code
fitted.model.bcfa.jags <- bcfa(
  model = cfa.model.syntax.lavaan,
  #dp=measurement.model.priors,
  n.chains=1,
  #burnin = n.burnin,
  burnin = 1,   # Think this is warmup in stan
  #adapt=n.warmup,
  sample=100,
  std.lv = TRUE,			# identify the model by fixing factor variance to 1
  #std.ov = TRUE,      # to get standardized solution
  int.ov.free = FALSE,
  #estimator = "ML",
  #likelihood = "wishart", # normal based  on dividing by N, wishart is N-1
  #se = "standard"
  target = "jags",
  mcmcfile=TRUE,
  save.lvs=FALSE,
  test="none",   # turn off computing fit functions
  #bcontrol=list(cores=n.chains),
  #mcmcextra = list(data='parvec=c(.5, .5, .5, .5, .5, .5, .5, .5)'),
  #mcmcextra = list(data=list(lambda=lambda.values)),
  mcmcextra = list(data=list(parvec=parvec.values)),
  inits = "jags",
  data = indicators
)

lambda.values <- array(NA, c(4,1,1))
lambda.values[,,1] <- .5
parvec.values=rep(.5, 8)





# Prepare to fit the structural model ------

# * Declare some MCMC parameters that hold over iterations from measurement model -----

# * * Initial values for structural model

if(model.has.outcomes){

  beta.o.inits <- NULL
  inv.psi.y.inits <- NULL
  #k=1

  for(k in 1:K){
    # Get initial values from regression of outcome on mean of indicators
    temp.y <- rowMeans(outcomes[k])
    temp.x <- rowMeans(indicators)

    temp.regression <- lm(temp.y ~ temp.x)

    beta.o.inits <- c(beta.o.inits, temp.regression$coefficients["temp.x"])
    inv.psi.y.inits <- c(inv.psi.y.inits, sigma(temp.regression)^2)

  } # closes loop over K

}

if(model.has.covariates){

  beta.c.inits <- NULL
  #h=1
  for(h in 1:H){
    # Get initial values from regression of mean of indicators on covariates
    temp.w <- rowMeans(covariates[h])
    temp.x <- rowMeans(indicators)

    temp.regression <- lm(temp.x ~ temp.w)
    summary(temp.regression)

    beta.c.inits <- c(beta.c.inits, temp.regression$coefficients["temp.w"])

  } # closes loop over H

}

# Define the initial values based on whether there are covariates or outcomes
if(model.has.outcomes & !model.has.covariates) {
  inits1 <- list(
    beta.o=beta.o.inits,
    inv.psi.y=inv.psi.y.inits
  )
}

if(!model.has.outcomes & model.has.covariates) {
  inits1 <- list(
    beta.c=beta.c.inits
  )
}

if(model.has.outcomes & model.has.covariates) {
  inits1 <- list(
    beta.c=beta.c.inits,
    beta.o=beta.o.inits,
    inv.psi.y=inv.psi.y.inits
  )
}














# Conduct MUPPET modeling with an outcome ------
MUPPET.CFA.outcome <- MUPPET.CFA.function(
  data=raw.data,
  indicators.names=c("PI", "FI", "AD", "FC"),
  covariates.names=NULL,
  outcomes.names="IGC",
  measurement.model.priors=dpriors(
    #nu="normal(0,3)",        # measured variable intercepts in N(mean, sd)
    lambda="normal(0,10)",    # measured variable loadings in N(mean, sd)
    #theta="gamma(.01,.01)[prec]"     # measured variable error precision
    theta="gamma(1,1)[prec]"     # measured variable error precision
    #psi="gamma(5,10)[prec]"        # latent variable precision
  ),
  combined.model.priors=dpriors(
    #itheta="dgamma(1,1)[prec]",     # measured variable error precision for indicators
    beta="dnorm(0,.0001)",     # structural coefficient in N(mean, precision)
    ipsi="dgamma(1,1)[prec]",     # variance terms for structural model (including error variance for outcomes)
    target="jags"
  ),
  n.chains = 2,
  n.warmup = 500,
  n.burnin = 0,
  n.iters.per.chain.after.warmup.and.burnin = 10,
  obtain.standardized.cfa=TRUE,
  obtain.standardized.combined.model = TRUE,
  beta.o.prior.mean = 0,
  beta.o.prior.var = 10000,
  psi.y.prior.alpha = 1,
  psi.y.prior.beta = 1,
  beta.c.prior.mean = 0,
  beta.c.prior.var = 10000,
  n.iters.per.chain.total.structural = 51,
  save.summary.stats.from.MUPPET=TRUE,
  save.draws.from.MUPPET=TRUE,
  model.check=FALSE,
  save.post.pred.data=FALSE
)












object <- fitted.model.bsem.jags

draws.from.bsem.object <- make_mcmc(fitted.model.bsem.jags@external$mcmcout)
draws.from.bsem.object <- do.call("rbind", draws.from.bsem.object)
n.draws.from.bsem <- nrow(draws.from.bsem.object)

replaced.draws <- rbind(
  draws.from.bsem.object,
  matrix(NA, nrow=nrow(draws.from.MUPPET.model)-n.draws.from.bsem, ncol=ncol(draws.from.bsem.object))
)


which.col=9
for(which.col in 1:ncol(replaced.draws)){
  column.name <- colnames(replaced.draws)[which.col]

  # check if the column name in the bsem object is also one for which there are draws from the combined model from jags
  # and if so, replace it
  # Otherwise, extend down the values of that column based on what's in the first row
  if(any(column.name == colnames(draws.from.MUPPET.model))){
    replaced.draws[ , column.name] <- draws.from.MUPPET.model[ , column.name]
  }
  if(!any(column.name == colnames(draws.from.MUPPET.model))){
    replaced.draws[ , column.name] <- replaced.draws[1, column.name]
  }
  #which(column.name == colnames(measurement.model.draws.as.data.frame))
}



standardizedPosterior(fitted.model.bsem.jags)

# Name the object with *most* of the correct bsem() object elements
object <- fitted.model.bsem.jags

#draws <- make_mcmc(object@external$mcmcout)
#draws <- do.call("rbind", draws)

# Used the replaced draws for summarizing
draws <- replaced.draws

# Rest of this is from standardizedPosterior() function
tf <- object
tmp <- fill_params(draws[1, ], object@Model, object@ParTable)

tf@Model <- tmp
tf@ParTable$est[tf@ParTable$free > 0] <- lav_model_get_parameters(tmp)

dots <- list(type="std.all", cov.std="TRUE")
tmp2 <- do.call("standardizedSolution", c(list(object = tf),
                                          dots))

#tmp2 <- do.call("standardizedSolution", c(list(object = tf)))

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



# MUPPET.CFA.outcome$summary.statistics.MUPPET
# MUPPET.CFA.outcome$ICDM.p.values

# Conduct MUPPET modeling with a covariate ------
MUPPET.CFA.covariate <- MUPPET.CFA.function(
  data=raw.data,
  indicators.names=c("PI", "FI", "AD", "FC"),
  covariates.names="age_in_years",
  outcomes.names=NULL,
  measurement.model.priors=dpriors(
    #nu="normal(0,3)",        # measured variable intercepts in N(mean, sd)
    lambda="normal(0,10)",    # measured variable loadings in N(mean, sd)
    #theta="gamma(.01,.01)[prec]"     # measured variable error precision
    theta="gamma(1,1)[prec]"     # measured variable error precision
    #psi="gamma(5,10)[prec]"        # latent variable precision
  ),
  n.chains = 2,
  n.warmup = 500,
  n.burnin = 0,
  n.iters.per.chain.after.warmup.and.burnin = 10000,
  obtain.standardized.cfa=TRUE,
  beta.o.prior.mean = 0,
  beta.o.prior.var = 10000,
  psi.y.prior.alpha = 1,
  psi.y.prior.beta = 1,
  beta.c.prior.mean = 0,
  beta.c.prior.var = 10000,
  n.iters.per.chain.total.structural = 51,
  save.summary.stats.from.MUPPET=TRUE,
  save.draws.from.MUPPET=TRUE,
  model.check=TRUE,
  save.post.pred.data=TRUE
)

# MUPPET.CFA.covariate$summary.statistics.MUPPET
# MUPPET.CFA.covariate$summary.statistics.PPMC

# MUPPET.CFA.covariate$LR.p.value
# MUPPET.CFA.covariate$resid.cor.p.values
# MUPPET.CFA.covariate$ICDM.p.values


