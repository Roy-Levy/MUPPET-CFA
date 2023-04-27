# Governing R code to run analyses for 
# Measurement and Uncertainty Preserving Parametric (MUPPET) CFA on
# A simulated example



# 1/25/2022
# This is my file, with folder structures and versions of the files for my development 
# This is not to be shared or made public. 
# The files to be shared are in the folder are in 
# E:\MUPPET CFA with Covariates and Outcomes\Simulated Example\Public
# Go there and run the file 'MUPPET CFA Simulated Example.R'


rm(list=ls())



# Define the folders ------

#   Drive letter
#   main analysis folder
#   R Code folder and functions subfolder
#   Conditions folder

#main.folder <- getwd()

library(rprojroot)
file.name.in.root.directory <- "MUPPET CFA Simulated Example.Rproj"
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
library(R2jags)
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
file.name <- "gen.data.1.dat"
raw.data <- as.data.frame(read.table(paste0(data.folder, file.name), header = TRUE))


# Conduct MUPPET modeling with an outcome ------
MUPPET.CFA.outcome <- MUPPET.CFA.function(
  data=raw.data,
  indicators.names=c("x1", "x2", "x3", "x4"),
  covariates.names=NULL,
  outcomes.names="y1",
  measurement.model.priors=dpriors(
    #nu="normal(0,3)",        # measured variable intercepts in N(mean, sd)
    lambda="normal(0,10)",    # measured variable loadings in N(mean, sd)
    #theta="gamma(.01,.01)[prec]"     # measured variable error precision
    theta="gamma(1,.8)[prec]"     # measured variable error precision
    #psi="gamma(5,10)[prec]"        # latent variable precision
  ),
  n.chains = 2,
  n.warmup = 500,
  n.burnin = 0,
  n.iters.per.chain.after.warmup.and.burnin = 1000,
  obtain.standardized.cfa=TRUE,
  beta.o.prior.mean = 0,
  beta.o.prior.var = 100,
  psi.y.prior.alpha = 1,
  psi.y.prior.beta = 1,
  beta.c.prior.mean = 0,
  beta.c.prior.var = 100,
  n.iters.per.chain.total.structural = 51,
  save.summary.stats.from.MUPPET=TRUE,
  save.draws.from.MUPPET=TRUE,
  model.check=TRUE,
  save.post.pred.data=TRUE
)


MUPPET.CFA.outcome$summary.statistics.MUPPET
MUPPET.CFA.outcome$ICDM.p.values

# Conduct MUPPET modeling with a covariate ------
MUPPET.CFA.covariate <- MUPPET.CFA.function(
  data=raw.data,
  indicators.names=c("x1", "x2", "x3", "x4"),
  covariates.names="w1",
  outcomes.names=NULL,
  measurement.model.priors=dpriors(
    #nu="normal(0,3)",        # measured variable intercepts in N(mean, sd)
    lambda="normal(0,10)",    # measured variable loadings in N(mean, sd)
    #theta="gamma(.01,.01)[prec]"     # measured variable error precision
    theta="gamma(1,.8)[prec]"     # measured variable error precision
    #psi="gamma(5,10)[prec]"        # latent variable precision
  ),
  n.chains = 2,
  n.warmup = 500,
  n.burnin = 0,
  n.iters.per.chain.after.warmup.and.burnin = 10000,
  obtain.standardized.cfa=TRUE,
  beta.o.prior.mean = 0,
  beta.o.prior.var = 100,
  psi.y.prior.alpha = 1,
  psi.y.prior.beta = 1,
  beta.c.prior.mean = 0,
  beta.c.prior.var = 100,
  n.iters.per.chain.total.structural = 51,
  save.summary.stats.from.MUPPET=TRUE,
  save.draws.from.MUPPET=TRUE,
  model.check=TRUE,
  save.post.pred.data=TRUE
)

MUPPET.CFA.covariate$summary.statistics.MUPPET
MUPPET.CFA.covariate$resid.cor.p.values
MUPPET.CFA.covariate$ICDM.p.values




