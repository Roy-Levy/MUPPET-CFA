# Governing R code to run analyses for
# Measurement and Uncertainty Preserving Parametric (MUPPET) CFA on
# The IIS Example
# Using the functions


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


# Conduct MUPPET modeling ------
IIS.measurement.model.blavaan.syntax <- '

  # measurement model
  # latent variable definitions
  f1 =~ NA*PI + NA*FI + NA*AD + NA*FC

  # mean structure (intercepts)
# 	PI ~ 0*1
# 	FI ~ 0*1
#   AD ~ 0*1
#   FC ~ 0*1

'


IIS.combined.model.with.outcome.blavaan.syntax <- '

  # measurement model
  # latent variable definitions
  f1 =~ NA*PI + NA*FI + NA*AD + NA*FC

  # mean structure (intercepts)
# 	PI ~ 0*1
# 	FI ~ 0*1
#   AD ~ 0*1
#   FC ~ 0*1

  # structural model
  IGC ~ f1

  # mean structure (intercepts)
  # IGC ~ 0*1
'

IIS.combined.model.with.covariate.blavaan.syntax <- '

  # measurement model
  # latent variable definitions
  f1 =~ NA*PI + NA*FI + NA*AD + NA*FC

  # mean structure (intercepts)
# 	PI ~ 0*1
# 	FI ~ 0*1
#   AD ~ 0*1
#   FC ~ 0*1

  # structural model
  f1 ~ age_in_years

  # mean structure (intercepts)
  f1 ~ 1
'



MUPPET.CFA.outcome <- MUPPET.CFA.function(
    data=raw.data,
    center.the.data=FALSE,
    #indicators.names=c("PI", "FI", "AD", "FC"),
    #covariates.names=NULL,
    #outcomes.names="IGC",
    measurement.model.blavaan.syntax=IIS.measurement.model.blavaan.syntax,
    combined.model.blavaan.syntax=IIS.combined.model.with.outcome.blavaan.syntax,
    combined.model.has.endogenous.latent=FALSE,
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
    n.iters.per.chain.after.warmup.and.burnin = 2500,
    obtain.standardized.measurement.model=TRUE,
    obtain.standardized.combined.model = TRUE,
    n.iters.per.chain.total.structural = 51,
    save.summary.stats.from.MUPPET=TRUE,
    save.draws.from.MUPPET=TRUE,
    model.check=FALSE,
    save.post.pred.data=FALSE
)

# MUPPET.CFA.outcome$summary.statistics.MUPPET




MUPPET.CFA.covariate<- MUPPET.CFA.function(
    data=raw.data,
    center.the.data=FALSE,
    #indicators.names=c("PI", "FI", "AD", "FC"),
    #covariates.names=NULL,
    #outcomes.names="IGC",
    measurement.model.blavaan.syntax=IIS.measurement.model.blavaan.syntax,
    combined.model.blavaan.syntax=IIS.combined.model.with.covariate.blavaan.syntax,
    combined.model.has.endogenous.latent=TRUE,
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
    n.iters.per.chain.after.warmup.and.burnin = 2500,
    obtain.standardized.measurement.model=TRUE,
    obtain.standardized.combined.model = TRUE,
    n.iters.per.chain.total.structural = 51,
    save.summary.stats.from.MUPPET=TRUE,
    save.draws.from.MUPPET=TRUE,
    model.check=FALSE,
    save.post.pred.data=FALSE
)


