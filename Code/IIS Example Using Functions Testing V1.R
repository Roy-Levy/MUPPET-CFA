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

#main.folder <- getwd()

library(rprojroot)
file.name.in.root.directory <- "MUPPET CFA IIS (Testing).Rproj"
main.folder <- find_root(has_file(file.name.in.root.directory))

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

file.name <- "MUPPET CFA function (mine).R"
source(paste0(functions.folder, file.name))


# Read in the data ------
file.name <- "IIS.dat"
raw.data <- as.data.frame(read.table(paste0(data.folder, file.name), header = TRUE))



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
n.iters.per.chain.after.warmup.and.burnin = 10000
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


# Frequentist joint modeling with a covariate ------

# * Define the model in the lavaan structure ------

lavaan.model <- '
  
  # measurement model
  # latent variable definitions
  f1 =~ NA*PI + NA*FI + NA*AD + NA*FC

	# variances
	f1 ~~ disturb_var*f1
	PI ~~ PI
	FI ~~ FI
	AD ~~ AD
	FC ~~ FC

  # means
  #f1 ~ 0*1
	#x1 ~ 0*1
	#x2 ~ 0*1
  #x3 ~ 0*1
  #x4 ~ 0*1


  # structural model
  f1 ~ beta*age_in_years

	# variances
  age_in_years ~~ covariate_var*age_in_years
	
  # means
  #w1 ~ 0*1
  
  ## MODEL CONSTRAINTS:
  disturb_var == 1 - (beta^2)*covariate_var
'


# * Fit the model ------

model.results <- sem(
  model = lavaan.model,
  data = raw.data,
  meanstructure = FALSE,
  orthogonal = FALSE,
  ##std.lv = TRUE,			# identify the model by fixing factor variance to 1
  std.lv = FALSE,			
  std.ov = FALSE,      # to get standardized solution
  estimator = "ML",
  likelihood = "wishart", # normal based  on dividing by N, wishart is N-1
  se = "standard"
  
)

# * Store the results -------------

# summary(model.results, standardized=FALSE)
# summary(model.results, standardized=TRUE)

# Obtain the standardized solution
# using "std.lv" forces the latent variable to have variance of 1
standardized.parameter.estimates <- standardizedSolution(
  model.results,
  #type="std.all",
  type="std.lv",
  se=TRUE,
  zstat=TRUE,
  pvalue=TRUE,
  ci=TRUE
) 

standardized.parameter.estimates <- standardizedSolution(
  model.results,
  type="std.all",
  #type="std.lv",
  se=TRUE,
  zstat=TRUE,
  pvalue=TRUE,
  ci=TRUE
) 



# Bayesian joint modeling with an outcome ------
# * Define a model in the lavaan structure ------


lavaan.model <- '
  
  # measurement model
  # latent variable definitions
  f1 =~ NA*PI + NA*FI + NA*AD + NA*FC

	# variances
	f1 ~~ factor_var*f1
	#f1 ~~ 1*f1
	PI ~~ PI
	FI ~~ FI
	AD ~~ AD
	FC ~~ FC

  # means
  f1 ~ 0*1
	PI ~ 0*1
	FI ~ 0*1
  AD ~ 0*1
  FC ~ 0*1


  # structural model
  IGC ~ f1

	# variances
  IGC ~~ IGC
	
  # means
  IGC ~ 0*1
  
  ## MODEL CONSTRAINTS:
  #factor_var == 1 
'



# * Define the prior distributions for families of parameters ------
if(1==1){
  
  model.priors.stan <- dpriors(
    #nu="normal(0,3)",        # measured variable intercepts in N(mean, sd)
    lambda="normal(0,10)",    # measured variable loadings in N(mean, sd)
    #theta="gamma(.01,.01)[prec]"     # measured variable error precision
    beta="normal(0,10)",    # measured variable loadings in N(mean, sd)
    theta="gamma(1,1)[prec]"     # measured variable error precision
    #psi="gamma(5,10)[prec]"        # latent variable precision
  )
  
}


# * Choose features of MCMC --------
#   the number of chains
#   the number of iterations to warmup
#	  the total number of iterations 
n.chains = 2
n.warmup = 500
n.burnin = 0
n.iters.per.chain.after.warmup.and.burnin = 10000
n.iters.total.per.chain = n.warmup+n.burnin+n.iters.per.chain.after.warmup.and.burnin

# * Fit the model ------

fitted.model.bsem <- bsem(
  model = lavaan.model,
  dp=model.priors.stan,
  n.chains=n.chains,
  #burnin = n.burnin,
  burnin = n.warmup,   # Think this is warmup in stan
  #adapt=n.warmup,
  sample=n.iters.per.chain.after.warmup.and.burnin,
  std.lv = TRUE,			# identify the model by fixing factor variance to 1
  #std.ov = TRUE,      # to get standardized solution
  #int.ov.free = FALSE,
  #estimator = "ML",
  #likelihood = "wishart", # normal based  on dividing by N, wishart is N-1
  #se = "standard"
  mcmcfile=FALSE,
  save.lvs=FALSE,
  test="none",   # turn off computing fit functions
  bcontrol=list(cores=n.chains),
  data = scale(raw.data, center = TRUE, scale = FALSE)
)


summary(fitted.model.bsem)


# * Obtain the standardized solution -----
standardized.solution <- standardizedPosterior(fitted.model.bsem)
# summary(standardized.solution)

# * Obtain summary statistics of the standardized solution -----
summary.statistics.standardized.solution <- MCMCsummary(
  standardized.solution,
  #params = parameters.to.summarize,
  HPD=TRUE,
  Rhat=FALSE, # won't compute because had to combine to one chain
  n.eff=FALSE, # won't compute because had to combine to one chain
  round=8, 
  func=median, 
  func_name = "median"
)

# * Save it as a new object  ------------
joint.model.with.outcome <- fitted.model.bsem


# Bayesian joint modeling with a covariate ------
# * Define a model in the lavaan structure ------


lavaan.model <- '
  
  # measurement model
  # latent variable definitions
  f1 =~ NA*PI + NA*FI + NA*AD + NA*FC

	# variances
	f1 ~~ factor_var*f1
	#f1 ~~ 1*f1
	PI ~~ PI
	FI ~~ FI
	AD ~~ AD
	FC ~~ FC

  # means
  #f1 ~ 0*1
	PI ~ 0*1
	FI ~ 0*1
  AD ~ 0*1
  FC ~ 0*1


  # structural model
  f1 ~ age_in_years

	# variances
  
	
  # means
  f1 ~ 0*1
  
  ## MODEL CONSTRAINTS:
  #factor_var == 1 
'


# * Define the prior distributions for families of parameters ------
if(1==1){
  
  model.priors.stan <- dpriors(
    #nu="normal(0,3)",        # measured variable intercepts in N(mean, sd)
    lambda="normal(0,10)",    # measured variable loadings in N(mean, sd)
    #theta="gamma(.01,.01)[prec]"     # measured variable error precision
    beta="normal(0,10)",    # measured variable loadings in N(mean, sd)
    theta="gamma(1,1)[prec]"     # measured variable error precision
    #psi="gamma(5,10)[prec]"        # latent variable precision
  )
  
}


# * Choose features of MCMC --------
#   the number of chains
#   the number of iterations to warmup
#	  the total number of iterations 
n.chains = 2
n.warmup = 500
n.burnin = 0
n.iters.per.chain.after.warmup.and.burnin = 10000
n.iters.total.per.chain = n.warmup+n.burnin+n.iters.per.chain.after.warmup.and.burnin

# * Fit the model ------

fitted.model.bsem <- bsem(
  model = lavaan.model,
  dp=model.priors.stan,
  n.chains=n.chains,
  #burnin = n.burnin,
  burnin = n.warmup,   # Think this is warmup in stan
  #adapt=n.warmup,
  sample=n.iters.per.chain.after.warmup.and.burnin,
  std.lv = TRUE,			# identify the model by fixing factor variance to 1
  #std.ov = TRUE,      # to get standardized solution
  #int.ov.free = FALSE,
  #estimator = "ML",
  #likelihood = "wishart", # normal based  on dividing by N, wishart is N-1
  #se = "standard"
  mcmcfile=FALSE,
  save.lvs=FALSE,
  test="none",   # turn off computing fit functions
  bcontrol=list(cores=n.chains),
  data = scale(raw.data, center = TRUE, scale = FALSE)
)


summary(fitted.model.bsem)


# * Obtain the standardized solution -----
standardized.solution <- standardizedPosterior(fitted.model.bsem)
# summary(standardized.solution)

# * Obtain summary statistics of the standardized solution -----
summary.statistics.standardized.solution <- MCMCsummary(
  standardized.solution,
  #params = parameters.to.summarize,
  HPD=TRUE,
  Rhat=FALSE, # won't compute because had to combine to one chain
  n.eff=FALSE, # won't compute because had to combine to one chain
  round=8, 
  func=median, 
  func_name = "median"
)

# * Save it as a new object  ------------
joint.model.with.covariate <- fitted.model.bsem


# Save the workspace ------
file.name <- "IIS Example Workspace.RData"
save.image(file=file.name)

