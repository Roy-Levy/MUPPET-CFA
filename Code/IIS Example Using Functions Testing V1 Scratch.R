# * Code for the the structural model -----
# * * Here, hard-coded for testing blavaan -----

combined.model.syntax.lavaan <- '
f1 =~ NA*x1 + NA*x2 + NA*x3 + NA*x4

x1 ~ 0*1
x2 ~ 0*1
x3 ~ 0*1
x4 ~ 0*1

y1 ~ f1
y1 ~ 0*1

'

fitted.model.bsem.jags <- bsem(
  model = combined.model.syntax.lavaan,
  #dp=measurement.model.priors,
  n.chains=2,
  #burnin = n.burnin,
  burnin = 1,   # Think this is warmup in stan
  #adapt=n.warmup,
  sample=2,
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
  #mcmcextra = list(data=list(parvec=parvec.values)),
  inits = "jags",
  data = cbind(indicators, outcomes)
)


# Extract the parameter table from blavaan for the CFA model using stan
combined.partable.jags <- as.data.frame(fitted.model.bsem.jags@ParTable)

combined.partable.jags <- rename(combined.partable.jags, parameter.name.jags = pxnames)


# cfa.partable.stan
# 
partable.join.mm.and.sm <- left_join(combined.partable.jags, cfa.partable.stan.estimated.parameters, by='parameter.name.jags')

partable.join.mm.and.sm <- (partable.join.mm.and.sm %>% filter(free.x >0))

partable.join.mm.and.sm <- rename(partable.join.mm.and.sm, parameter.number.cm = free.x)
partable.join.mm.and.sm <- rename(partable.join.mm.and.sm, parameter.number.mm = stanpnum)

n.estimated.parameters.cm <- nrow(partable.join.mm.and.sm)
n.estimated.parameters.mm <- max(select(partable.join.mm.and.sm, parameter.number.mm), na.rm=TRUE)


parvec.values <- rep(NA, n.estimated.parameters.cm)
for(which.mm.param in 1:n.estimated.parameters.mm){
  place.in.parvec.values <- partable.join.mm.and.sm %>% 
    filter(parameter.number.mm == which.mm.param)  %>%
    select(parameter.number.cm)
  
  parvec.values[as.integer(place.in.parvec.values)] = .5
}



fitted.model.bsem.jags <- bsem(
  model = combined.model.syntax.lavaan,
  #dp=measurement.model.priors,
  n.chains=2,
  #burnin = n.burnin,
  burnin = 1,   # Think this is warmup in stan
  #adapt=n.warmup,
  sample=4,
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
 # mcmcextra = list(data=list(parvec=parvec.values)),
 # mcmcextra = list(data=list(parvec=rep(NA,10))),
 # mcmcextra = list(data=list(parvec=c(1,rep(NA,9)))),
  inits = "jags",
  data = cbind(indicators, outcomes)
)


# Load the object called 'jagtrans' with information for jags
load(
  paste0(getwd(), "/lavExport/semjags.rda")
)


library(R2jags)

# jags.data <- list("x"=x, "J"=J, "n"=n, "y"=y, "lambda"=lambda, "inv.psi.x"=inv.psi.x, "phi"=phi)
parvec.values <- rep(NA, 10)
parvec.values <- c(1, rep(NA, 9))

jags.data <- jagtrans$data


jags.data[[length(jags.data)+1]] <- parvec.values
names(jags.data)[[length(jags.data)]] <- "parvec"

# * Fit the model ------

# Load modules for JAGS
load.module("glm")

model.file.name <- "sem.jag"
model.file.name <- paste0(getwd(), "/lavExport/sem.jag")

# Initialize the model
jags.model.initialized <- jags.model(file=model.file.name,
                                     data=jags.data,
                                     #data=jags.data.no.factor.variance,
                                     n.chains=n.chains,
                                     #inits=inits
                                     #inits=inits1
)

# Now obtain the distribution
jags.model.fitted <- coda.samples(
  jags.model.initialized,
  variable.names=c("parvec","lambda"),
  #variable.names=entities.to.monitor,
  #variable.names=entities.to.monitor.factor.variance,
  #n.iter=n.iters.total.per.chain,
  #n.iter=n.iters.per.chain.total.structural,
  n.iter=10,
  #progress.bar="gui"
  progress.bar="none"
)



fitted.model.bsem.jags <- modified.bsem.fuction(
  model = combined.model.syntax.lavaan,
  #dp=measurement.model.priors,
  n.chains=2,
  #burnin = n.burnin,
  burnin = 1,   # Think this is warmup in stan
  #adapt=n.warmup,
  sample=10,
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
  # mcmcextra = list(data=list(parvec=parvec.values)),
  # mcmcextra = list(data=list(parvec=rep(NA,10))),
  # mcmcextra = list(data=list(parvec=c(1,rep(NA,9)))),
  inits = "jags",
  data = cbind(indicators, outcomes)
)

# summary(fitted.model.bsem.jags)

fitted.model.jags <- blavInspect(fitted.model.bsem.jags, what="mcobj")
(draws.as.mcmc.list <- coda::as.mcmc.list(fitted.model.jags))

fitted.model.jags$data
fitted.model.jags$model



# convert draws to mcmc.list
draws.as.mcmc.list <- MCMCchains(fitted.model.jags, 
                                 mcmc.list = TRUE)


# Define the draws to analyze
draws.to.analyze <- draws.as.mcmc.list



head((as.matrix(draws.to.analyze)))