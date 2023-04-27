# Governing R code to run analyses for 
# Measurement and Uncertainty Preserving Parametric (MUPPET) CFA on
# The IIS Example
# Using the functions


# 1/24/2022


rm(list=ls())



# Define the folders ------

#   Drive letter
#   main analysis folder
#   R Code folder and functions subfolder
#   Conditions folder

#main.folder <- getwd()

library(rprojroot)
file.name.in.root.directory <- "MUPPET CFA IIS.Rproj"
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

