# Function to conduct MUPPET CFA analyses

MUPPET.CFA.function <- function(
    data,
    center.the.data=FALSE,
    indicators.names,
    covariates.names=NULL,
    outcomes.names=NULL,
    measurement.model.blavaan.syntax=NULL,
    combined.model.blavaan.syntax=NULL,
    measurement.model.priors=blavaan::dpriors(),
    combined.model.priors=blavaan::dpriors(target="jags"),
    combined.model.has.endogenous.latent=FALSE,
    n.chains = 2,
    n.warmup = 500,
    n.burnin = 0,
    n.iters.per.chain.after.warmup.and.burnin = 10000,
    obtain.standardized.measurement.model=TRUE,
    obtain.standardized.combined.model=TRUE,
    # beta.o.prior.mean = 0, deprecated, now use combined.model.priors
    # beta.o.prior.var = 10000, deprecated, now use combined.model.priors
    # psi.y.prior.alpha = 1, deprecated, now use combined.model.priors
    # psi.y.prior.beta = 1, deprecated, now use combined.model.priors
    beta.c.prior.mean = 0,
    beta.c.prior.var = 10000,
    n.iters.per.chain.total.structural = 51,
    save.summary.stats.from.MUPPET=TRUE,
    save.draws.from.MUPPET=FALSE,
    model.check=FALSE,
    save.post.pred.data=FALSE
  ){

    # Define the data to analyze, possible to be replaced by centered data ------
    data.to.analyze <- data

    # Center the data if desired------
    if(center.the.data==TRUE){
      data.centered <- as.data.frame(scale(data.to.analyze, scale=FALSE))
      data.to.analyze <- data.centered
    }

    # STAGE 2: MEASUREMENT MODEL MODEL ----


    # Define the number of iterations per chain ------
    n.iters.total.per.chain = n.warmup+n.burnin+n.iters.per.chain.after.warmup.and.burnin

    # Print out beginning of fitting this portion of the model
    print(paste0("Fitting the measurement model for ",
                 n.chains, " chains and ",
                 n.iters.total.per.chain, " total iterations per chain")
    )

    # Fit the measurement model ------

    fitted.model.bcfa <- bcfa(
      #model = cfa.model.syntax.lavaan,
      model = measurement.model.blavaan.syntax,
      dp=measurement.model.priors,
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
      #data = indicators
      #data = data.centered
      data = data.to.analyze
    )

    # summary(fitted.model.bcfa)


    # Obtain the standardized solution to the measurement model -----
    if(obtain.standardized.measurement.model) standardized.posterior.measurement.model <- standardizedPosterior(fitted.model.bcfa)


    # STAGE 2: COMBINED MODEL ----

    # Extract the stanfit object for the fitted measurement model -----
    fitted.model.stanfit <- blavInspect(fitted.model.bcfa, what="mcobj")

    # convert draws to mcmc.list
    draws.as.mcmc.list <- MCMCchains(fitted.model.stanfit,
                                     mcmc.list = TRUE)

    # Define the draws to analyze
    draws.to.analyze <- draws.as.mcmc.list

    # Extract the parameter table from blavaan for the CFA model using stan
    measurement.model.partable.stan <- as.data.frame(fitted.model.bcfa@ParTable)

    # Define the internal names of the estimated parameters in stan
    parameter.name.pxnames <- (measurement.model.partable.stan %>% filter(!is.na(pxnames)))$pxnames

    # # Convert the draws from the measurement model to a data frame ------
    measurement.model.draws.as.data.frame <- dplyr::select(
      as.data.frame(as.matrix(draws.to.analyze)),
      contains(parameter.name.pxnames)
    )

    # Get what we think blavaan will call the parameters in JAGS, from the stan results -----

    # Initialize a vector
    parameter.name.jags <- c()

    # Obtain just the parameters that get a parameter number in blavaan
    temp <- measurement.model.partable.stan %>% filter(!is.na(stanpnum))

    # For each parameter, get the name in jags
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

    # Append the jags names for the parameter to the parameter table
    measurement.model.partable.stan.estimated.parameters <- cbind(temp, parameter.name.jags)

    # Rename the columns for the draws by their jags names
    colnames(measurement.model.draws.as.data.frame) <- parameter.name.jags

    # Prepare to fit the structural model ------

    # * Declare if standardizing the latent variable in fitting the combined model  ------
    # Default: yes for models without endogenous latent variable, no for models with endogenous latent variable
    # Controlled by the argument 'combined.model.has.endogenous.latent'

      standardize.lv.when.fitting.combined.model = !combined.model.has.endogenous.latent

      # * Run bsem() for a few iterations, just to produce JAGS syntax and data structure ------

      fitted.model.bsem.jags <- bsem(
        # model = combined.model.syntax.lavaan,
        model = combined.model.blavaan.syntax,
        dp=combined.model.priors,
        n.chains=2,
        #burnin = n.burnin,
        burnin = 1,   # Think this is warmup in stan
        adapt=1,
        sample=2,

        #std.lv = TRUE,			# identify the model by fixing factor variance to 1
        std.lv = standardize.lv.when.fitting.combined.model,			# identify the model by fixing factor variance to 1

        #std.ov = TRUE,      # to get standardized solution
        #int.ov.free = FALSE,
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
        # data = cbind(indicators, outcomes)
        #data = data.centered
        data = data.to.analyze
      )

      # * Define the file name with the JAGS syntax ------
      model.file.name <- paste0(getwd(), "/lavExport/sem.jag")


      # * Extract the parameter table from blavaan for the combined model using jags -----
      combined.partable.jags <- as.data.frame(fitted.model.bsem.jags@ParTable)
      combined.partable.jags <- rename(combined.partable.jags, parameter.name.jags = pxnames)

      # * Create the parameter table for the combined model -----
      # Merge the parameter table from the measurement model and combined model
      partable.join.mm.and.cm <- left_join(combined.partable.jags, measurement.model.partable.stan.estimated.parameters, by='parameter.name.jags')

      # Only preserve the rows for paramters that are estimated
      partable.join.mm.and.cm <- (partable.join.mm.and.cm %>% filter(free.x >0))

      # Rename the columns for the parameter number in each model
      partable.join.mm.and.cm <- rename(partable.join.mm.and.cm, parameter.number.cm = free.x)
      partable.join.mm.and.cm <- rename(partable.join.mm.and.cm, parameter.number.mm = stanpnum)

      # * Define the number of parameters in the combined model and measurement model -----
      n.estimated.parameters.cm <- nrow(partable.join.mm.and.cm)
      n.estimated.parameters.mm <- max(select(partable.join.mm.and.cm, parameter.number.mm), na.rm=TRUE)

      # * Load the object called 'jagtrans' with information for jags -----
      load(
        paste0(getwd(), "/lavExport/semjags.rda")
      )

      # * Extract the data portion that was passed to JAGS, to be augmented later -----
      jags.data <- jagtrans$data


      # * Change values for draws from earlier stage to be passed to JAGS via parvec -----
      # This has to be done because JAGS uses precision parameterization, while results from earlier stages are variances
      # The values in the variance columns terms created here are therefore not variances, but precisions
      measurement.model.draws.as.data.frame.changed <- measurement.model.draws.as.data.frame
      colnames(measurement.model.draws.as.data.frame.changed)

      # * * Read in the text file with the JAGS code ----------
      jags.code.as.vector <- readLines(model.file.name)

      # * * Invert values for variance terms ----------

      # Go through each parameter
      #which.param=5
      for(which.param in 1:ncol(measurement.model.draws.as.data.frame.changed)){

        # Get the parameter name
        parameter.name <- colnames(measurement.model.draws.as.data.frame.changed)[which.param]

        # See if that parameter name is involved in inverting a value from parvec
        # Using first 2 characters from parameter name, to get rid of []
        # and seeing if
        #    the text 'parvec' is there and
        #    if there's a "-1" to indicate inverting

        text.to.search.for <- c(substr(parameter.name, 1, 2), "parvec", "-1")


        num.lines.with.text.to.search.for <- sum(
          rowSums(sapply(X = text.to.search.for, FUN = grepl, jags.code.as.vector))==length(text.to.search.for)
        )

        # If that parameter is involved, invert it
        if(num.lines.with.text.to.search.for > 0){
          measurement.model.draws.as.data.frame.changed[ , parameter.name] = 1/measurement.model.draws.as.data.frame.changed[ , parameter.name]
        }

      }

      # * Define entities to monitor ----------
      # entities.to.monitor <- c("beta.o", "beta.o.standardized.mod.imp.var.y", "psi.y", "lambda", "phi", "psi.x")
      # entities.to.monitor <- c("lambda", "theta", "beta", "psi")

      # Here, by the entities that were monitored in the earlier run of the combined model
      entities.to.monitor <- fitted.model.bsem.jags@external$mcmcout$monitor

      # Fit the structural portion of the model --------
      print(
        paste0("Fitting the structural model for ", nrow(measurement.model.draws.as.data.frame), " draws for measurement model")
      )


      # * Call the packages ------
      #library(parallel)
      #library(foreach)
      #library(doParallel)


      # * Clean up existing parallel environment ------
      unregister_dopar <- function() {
        env <- foreach:::.foreachGlobals
        rm(list=ls(name=env), pos=env)
      }

      unregister_dopar()


      # * Setup parallel environment ------

      # * * Detect the number of cores
      numCores <- detectCores()

      # * * Set the number of cores to use
      registerDoParallel(numCores)  # use multicore, set to the number of our cores


      # * Run analyses in parallel ------
      # which.iter=1
      #results.from.parallel <- foreach(which.iter=1:10,
      draws.from.MUPPET.model <- foreach(which.iter=1:nrow(measurement.model.draws.as.data.frame),
                                         #results.from.parallel <- foreach(which.iter=1:1000,
                                         .packages = c("dplyr", "R2jags"),
                                         .combine=rbind) %dopar% {


                                           # * * Create the vector of values for parameters to be fed to JAGS -----
                                           parvec.values.for.iter <- rep(NA, n.estimated.parameters.cm)
                                           for(which.mm.param in 1:n.estimated.parameters.mm){
                                             place.in.parvec.values <- partable.join.mm.and.cm %>%
                                               filter(parameter.number.mm == which.mm.param)  %>%
                                               select(parameter.number.cm)

                                             # parvec.values.for.iter[as.integer(place.in.parvec.values)] = measurement.model.draws.as.data.frame[which.iter, which.mm.param]
                                             # Changing to now use the *changed* values of the measurement model from above, which inverts the variance terms for passing to parvec
                                             # Because JAGS uses precision parameterization
                                             parvec.values.for.iter[as.integer(place.in.parvec.values)] = measurement.model.draws.as.data.frame.changed[which.iter, which.mm.param]

                                           }

                                           # * * Augment the data to be passed to jags with parvec values just selected -----
                                           jags.data.with.parvec <- jags.data
                                           jags.data.with.parvec[[length(jags.data)+1]] <- parvec.values.for.iter
                                           names(jags.data.with.parvec)[[length(jags.data.with.parvec)]] <- "parvec"

                                           #jags.data[[length(jags.data)+1]] <- parvec.values.for.iter
                                           #names(jags.data)[[length(jags.data)]] <- "parvec"


                                           # * * Fit the model ------

                                           # Load modules for JAGS
                                           load.module("glm")

                                           # Initialize the model
                                           jags.model.initialized <- jags.model(file=model.file.name,
                                                                                data=jags.data.with.parvec,
                                                                                #data=jags.data,
                                                                                #data=jags.data.no.factor.variance,
                                                                                n.chains=n.chains,
                                                                                #inits=inits
                                                                                #inits=inits1
                                           )



                                           # Now obtain the distribution
                                           jags.model.fitted <- coda.samples(
                                             jags.model.initialized,

                                             # variable.names=c(entities.to.monitor, "parvec"),
                                             variable.names=entities.to.monitor,


                                             #variable.names=entities.to.monitor.factor.variance,
                                             #n.iter=n.iters.total.per.chain,
                                             n.iter=n.iters.per.chain.total.structural,
                                             #progress.bar="gui"
                                             progress.bar="none"
                                           )




                                           # * * Store the value from the desired iteration ------
                                           result1 <- jags.model.fitted[n.iters.per.chain.total.structural, ][[1]]
                                           result1


                                         } # closes foreach


      # * Close parallel environment ------
      stopImplicitCluster()





    # STAGE 3: MODEL CHECKING ----
    if(model.check){

      # Fit the structural portion of the model --------
      print(
        paste0("Conducting model checking for ", nrow(draws.from.MUPPET.model), " draws")
      )

      # Call the packages -----
      #library(MASS)
      #library(lavaan)

      # Set up folder for this analysis ------
      PPMC.folder <- paste0(getwd(), "/PPMC/")
      dir.create(PPMC.folder)
      #setwd(PPMC.folder)

      # Define the draws to use from the model ------
      draws <- as.data.frame(draws.from.MUPPET.model)

      # Clean up existing parallel environment ------
      unregister_dopar <- function() {
        env <- foreach:::.foreachGlobals
        rm(list=ls(name=env), pos=env)
      }

      unregister_dopar()

      # Setup parallel environment ------

      # Detect the number of cores
      numCores <- detectCores()

      # Set the number of cores to use
      registerDoParallel(numCores)  # use multicore, set to the number of our cores

      # Conduct PPMC ------

      # which.iter = 1
      # which.iter = which.iter+1
      PPMC.results.from.parallel <- foreach(which.iter=1:nrow(draws),
                                            .export = c("LR.function", "resid.cor.matrix.function", "ICDM.function"),
                                            .packages = c("dplyr", "MASS", "rstatix"),
                                            .combine=rbind) %dopar% {


                                              PPMC.result = NULL

                                              # Set up numbers of phony variables for covariates and outcomes if they're not there
                                              if(!model.has.covariates) H=1
                                              if(!model.has.outcomes) K=1


                                              # * GENERATE MODEL-IMPLIED COVARIANCE MATRIX AND POSTERIOR PREDICTED DATA -----

                                              # * * Define the model-implied covariance matrix based representation -----
                                              # This is essentially treating the covariate as a single-indicator of a latent variable
                                              # This is essentially treating the outcome as a single-indicator of a latent variable

                                              # * * * The loading matrix -----
                                              # Has 1 for the covariates and outcomes (as a single indicator structure)
                                              # and value of 'loading' for indicator
                                              loading.matrix.combined <- NULL

                                              for(h in 1:H){
                                                loading.matrix.combined <- rbind(loading.matrix.combined, c(rep(1, H), rep(0,M), rep(0,K)))
                                              }

                                              # loadings for indicators
                                              # select draws based on colname having 'lambda'
                                              temp <- draws[
                                                which.iter,
                                                grep(
                                                  'lambda',
                                                  colnames(draws),
                                                  value=TRUE
                                                )
                                              ]

                                              for(j in 1:J){
                                                loading.matrix.combined <- rbind(loading.matrix.combined, c(rep(0, H), as.numeric(temp[j]), rep(0,K)))
                                              }

                                              for(k in 1:K){
                                                loading.matrix.combined <- rbind(loading.matrix.combined, c(rep(0, H), rep(0,M), rep(1,K)))
                                              }


                                              # * * * The structural coefficient matrix ------
                                              # Has 0 for the first row (representing the implicit latent variable for the covariate in a single indicator setup)
                                              # Other values extract structural coefficient values depending on whether the model has covariates or outcomes

                                              beta.matrix <- matrix(0, nrow=H+M+K, ncol=H+M+k)

                                              if(model.has.outcomes){
                                                # select draws based on colname having 'beta.o'
                                                temp <- draws[
                                                  which.iter,
                                                  grep(
                                                    'beta.o',
                                                    colnames(draws),
                                                    value=TRUE
                                                  )
                                                ]

                                                # Now get rid of the ones that are standardized
                                                temp <- temp[
                                                  ,
                                                  grep(
                                                    'standardized',
                                                    colnames(temp),
                                                    value=TRUE,
                                                    invert=TRUE
                                                  )
                                                ]

                                                beta.matrix[3, 2] <- as.numeric(temp)


                                              } # closes if outcomes present

                                              if(model.has.covariates){
                                                # select draws based on colname having 'beta.c'
                                                temp <- draws[
                                                  which.iter,
                                                  grep(
                                                    'beta.c',
                                                    colnames(draws),
                                                    value=TRUE
                                                  )
                                                ]

                                                # Now get rid of the ones that are standardized
                                                temp <- temp[
                                                  ,
                                                  grep(
                                                    'standardized',
                                                    colnames(temp),
                                                    value=TRUE,
                                                    invert=TRUE
                                                  )
                                                ]

                                                beta.matrix[2, 1] <- as.numeric(temp)


                                              } # closes if covariates present


                                              # * * * The error covariance matrix -----
                                              # Has 0 for the variance of the covariate (single indicator model)
                                              # Has 0 for the variance of the outcome (single indicator model)
                                              #
                                              error.matrix.combined <- matrix(0, nrow=H+J+K, ncol=H+J+K)

                                              # select draws based on colname having 'psi.x'
                                              temp <- draws[
                                                which.iter,
                                                grep(
                                                  'psi.x',
                                                  colnames(draws),
                                                  value=TRUE
                                                )
                                              ]

                                              for(j in 1:J){
                                                error.matrix.combined[j+H,j+H] <- as.numeric(temp[1,j])
                                              }


                                              # * * * The latent disturbance covariance matrix -----
                                              # Has the variance for the covariate (single indicator model)
                                              # and value of '1 -  squared structural.w' for the endogenous factor
                                              # and value of 'variances.y -  squared structural.y' for the outcome (single indicator model)

                                              # latent variables for covariates
                                              if(model.has.covariates) disturbance.vector.latent.for.covariates <- diag(var(covariates))
                                              if(!model.has.covariates) disturbance.vector.latent.for.covariates <- diag(1, H)

                                              # latent variables for indicators
                                              if(model.has.covariates){
                                                # select draws based on colname having 'psi.f'
                                                temp <- draws[
                                                  which.iter,
                                                  grep(
                                                    'psi.f',
                                                    colnames(draws),
                                                    value=TRUE
                                                  )
                                                ]

                                                disturbance.vector.latent.for.indicators <- as.numeric(temp)
                                              }
                                              if(!model.has.covariates){
                                                # select draws based on colname having 'phi'
                                                temp <- draws[
                                                  which.iter,
                                                  grep(
                                                    'phi',
                                                    colnames(draws),
                                                    value=TRUE
                                                  )
                                                ]

                                                disturbance.vector.latent.for.indicators <- as.numeric(temp)
                                              }

                                              # latent variables for outcomes
                                              if(model.has.outcomes){
                                                # select draws based on colname having 'psi.y'
                                                temp <- draws[
                                                  which.iter,
                                                  grep(
                                                    'psi.y',
                                                    colnames(draws),
                                                    value=TRUE
                                                  )
                                                ]

                                                disturbance.vector.latent.for.outcomes <- as.numeric(temp)

                                              }
                                              if(!model.has.outcomes) disturbance.vector.latent.for.outcomes <- diag(1, K)

                                              # Put the pieces together for the latent disturbance covariance matrix
                                              disturbance.matrix.combined <- diag(c(
                                                disturbance.vector.latent.for.covariates,
                                                disturbance.vector.latent.for.indicators,
                                                disturbance.vector.latent.for.outcomes
                                              ))

                                              # * * * Put the pieces together to form the covariance matrix of observables -----
                                              # See chapter by Bollen in Handbook by Boyle, equation 2.13
                                              # I think that's missing a transposed loading matrix at the end
                                              # See the first element in equation 2.22. That has the transposed loading matrix
                                              # The order is
                                              #   all the covariates
                                              #   all the indicators
                                              #   all the outcomes
                                              mod.imp.cov.matrix.observables <-
                                                loading.matrix.combined %*%
                                                solve((diag(1, nrow=ncol(beta.matrix)) - beta.matrix)) %*%
                                                disturbance.matrix.combined %*%
                                                t(solve((diag(1, nrow=ncol(beta.matrix)) - beta.matrix))) %*%
                                                t(loading.matrix.combined) +
                                                error.matrix.combined

                                              # Add column names
                                              column.names <- c(
                                                paste("w", seq(1:H), sep=""),
                                                paste("x", seq(1:J), sep=""),
                                                paste("y", seq(1:K), sep="")
                                              )

                                              colnames(mod.imp.cov.matrix.observables) <- column.names

                                              # * * * Select portion of model-implied covariance matrix based on whether there are covariates or outcomes ----
                                              if(model.has.covariates & !model.has.outcomes) mod.imp.cov.matrix.observables <- mod.imp.cov.matrix.observables[1:(H+J), 1:(H+J)]
                                              if(!model.has.covariates & model.has.outcomes) mod.imp.cov.matrix.observables <- mod.imp.cov.matrix.observables[(H+1):(H+J+K), (H+1):(H+J+K)]

                                              # * * Generate the data -----
                                              if(model.has.covariates & !model.has.outcomes) post.pred.data.raw <- mvrnorm(n=n, mu=rep(0, H+J), Sigma=mod.imp.cov.matrix.observables)
                                              if(!model.has.covariates & model.has.outcomes) post.pred.data.raw <- mvrnorm(n=n, mu=rep(0, J+K), Sigma=mod.imp.cov.matrix.observables)
                                              if(model.has.covariates & model.has.outcomes) post.pred.data.raw <- mvrnorm(n=n, mu=rep(0, H+J+K), Sigma=mod.imp.cov.matrix.observables)

                                              # * * * Center so that means are 0 -----
                                              post.pred.data.centered <- scale(post.pred.data.raw, center=TRUE, scale=FALSE)

                                              # * * * Use the centered data ------
                                              post.pred.data <- post.pred.data.centered

                                              if(model.has.covariates & !model.has.outcomes) colnames(post.pred.data) <- column.names[1:(H+J)]
                                              if(!model.has.covariates & model.has.outcomes) colnames(post.pred.data) <- column.names[(H+1):(H+J+K)]
                                              if(model.has.covariates & model.has.outcomes) colnames(post.pred.data) <- column.names[1:(H+J+K)]

                                              # * * Write out the data ------
                                              if(save.post.pred.data){
                                                data.file.name <- paste0("post.pred.data.iter.", which.iter, ".dat")
                                                write.table(
                                                  x=post.pred.data,
                                                  file=paste0(PPMC.folder, data.file.name),
                                                  row.names=FALSE,
                                                  col.names=TRUE
                                                )
                                              } # closes if saving posterior predictive data



                                              # * Define the observables to analyze -----
                                              observables.to.analyze <- indicators
                                              if(model.has.covariates) observables.to.analyze <- cbind(covariates, observables.to.analyze)
                                              if(model.has.outcomes) observables.to.analyze <- cbind(observables.to.analyze, outcomes)


                                              # * Likelihood Ratio Discrepancy Measure PPMC -----
                                              # * * Compute the realized value ----
                                              realized.LR.fit <- LR.function(
    data.cov.matrix=cov(observables.to.analyze),
    mod.imp.cov.matrix=mod.imp.cov.matrix.observables,
    n=n
                                              )

                                              # * * Compute the posterior predicted value ----
                                              postpred.LR.fit <- LR.function(
    data.cov.matrix=cov(post.pred.data),
    mod.imp.cov.matrix=mod.imp.cov.matrix.observables,
    n=n
                                              )

                                              # * * Bind the realized and posterior predicted values ----
                                              realized.and.post.pred.values <- c(realized.LR.fit, postpred.LR.fit)
                                              names(realized.and.post.pred.values) <- colnames(cbind(realized.LR.fit, postpred.LR.fit))

                                              # * * Store the result as part of possibly multiple results
                                              PPMC.result <- c(PPMC.result, realized.and.post.pred.values)


                                              # * Residual Correlation Discrepancy Measure PPMC -----

                                              # * * Compute the realized values ----
                                              realized.resid.cor.wide <- resid.cor.matrix.function(
    data.matrix = observables.to.analyze,
    mod.imp.cov.matrix = mod.imp.cov.matrix.observables,
    output="wide"
                                              )
                                              names(realized.resid.cor.wide) <- paste0("realized_resid_cor.", names(realized.resid.cor.wide))

                                              # * * Compute the posterior predicted values ----
                                              postpred.resid.cor.wide <- resid.cor.matrix.function(
    data.matrix = post.pred.data,
    mod.imp.cov.matrix = mod.imp.cov.matrix.observables,
    output="wide"
                                              )
                                              names(postpred.resid.cor.wide) <- paste0("postpred_resid_cor_", names(postpred.resid.cor.wide))

                                              # * * Bind the realized and posterior predicted values ----
                                              realized.and.post.pred.values <- c(realized.resid.cor.wide, postpred.resid.cor.wide)

                                              # * * Store the result as part of multiple results
                                              PPMC.result <- c(PPMC.result, realized.and.post.pred.values)



                                              # * ICDM Discrepancy Measure PPMC -----

                                              # * * Extract the loadings for indicators from this iteration ----
                                              # select draws based on colname having 'lambda'
                                              loadings.this.iter <- draws[
                                                which.iter,
                                                grep(
                                                  'lambda',
                                                  colnames(draws),
                                                  value=TRUE
                                                )
                                              ]

                                              # * * Define the external (non-indicator) variables  ----
                                              if(model.has.covariates & !model.has.outcomes) non.indicators.data.matrix = covariates
                                              if(!model.has.covariates & model.has.outcomes) non.indicators.data.matrix = outcomes
                                              if(model.has.covariates & model.has.outcomes) non.indicators.data.matrix = cbind(covariates, outcomes)


                                              # * * Compute the realized values ----
                                              realized.ICDM.wide <- ICDM.function(
    indicators.data.matrix = indicators,
    indicators.loadings = loadings.this.iter,
    non.indicators.data.matrix = non.indicators.data.matrix
                                              )
                                              colnames(realized.ICDM.wide) <- paste0("realized_", colnames(realized.ICDM.wide))

                                              # * * Compute the posterior predicted values ----
                                              # * * * Define the posterior predicted data matrix for just the indicators and externals ----
                                              # select columns based on colname having 'x'
                                              post.pred.data.indicators <- post.pred.data[
                                                ,
                                                grep(
                                                  'x',
                                                  colnames(post.pred.data),
                                                  value=TRUE
                                                )
                                              ]

                                              # select columns based on colname not having 'x'
                                              post.pred.data.non.indicators <- dplyr::select(as.data.frame(post.pred.data), !contains("x"))

                                              # * * Now compute the posterior predicted values ----
                                              postpred.ICDM.wide <- ICDM.function(
    indicators.data.matrix = post.pred.data.indicators,
    indicators.loadings = loadings.this.iter,
    non.indicators.data.matrix = post.pred.data.non.indicators
                                              )
                                              colnames(postpred.ICDM.wide) <- paste0("postpred_", colnames(postpred.ICDM.wide))


                                              # * * Bind the realized and posterior predicted values ----
                                              realized.and.post.pred.values <- cbind(realized.ICDM.wide, postpred.ICDM.wide)
                                              names <- colnames(realized.and.post.pred.values)
                                              realized.and.post.pred.values <- as.vector(realized.and.post.pred.values)
                                              names(realized.and.post.pred.values) <- names

                                              # * * Store the result as part of multiple results
                                              # PPMC.result$ICDM.result <- realized.and.post.pred.values
                                              PPMC.result <- c(PPMC.result, realized.and.post.pred.values)


                                              return(PPMC.result)

                                            } # closes foreach

      # * Close parallel environment ------
      stopImplicitCluster()



      # * Compute p-values ------
      # * * LR ----
      # * * * Select PPMC results for discrepancy measure of interest
      temp <- PPMC.results.from.parallel[
        ,
        grep(
          'LR',
          colnames(PPMC.results.from.parallel),
          value=TRUE
        )
      ]

      p.value <- mean(temp[ ,grep(
        'postpred',
        colnames(temp),
        value=TRUE
      )
      ] >=
        temp[ ,grep(
          'realized',
          colnames(temp),
          value=TRUE
        )
        ]
      )

      LR.p.value <- p.value



      # * * Residual Correlation ----
      # * * * Select PPMC results for discrepancy measure of interest

      temp <- PPMC.results.from.parallel[
        ,
        grep(
          'resid',
          colnames(PPMC.results.from.parallel),
          value=TRUE
        )
      ]

      # * * * Isolate realized values
      temp.realized.values <- temp[
        ,
        grep("realized", colnames(temp), value=TRUE)
      ]

      # * * * Isolate postpred values
      temp.postpred.values <- temp[
        ,
        grep("postpred", colnames(temp), value=TRUE)
      ]

      # * * * Compute p-values
      resid.cor.p.values <- colMeans(temp.realized.values <= temp.postpred.values)




      # * * ICDM Discrepancy Measure ----
      # * * * Select PPMC results for discrepancy measure of interest
      temp <- PPMC.results.from.parallel[
        ,
        grep(
          'ICDM',
          colnames(PPMC.results.from.parallel),
          value=TRUE
        )
      ]

      # * * * Isolate realized values
      temp.realized.values <- temp[
        ,
        grep("realized", colnames(temp), value=TRUE)
      ]

      # * * * Isolate postpred values
      temp.postpred.values <- temp[
        ,
        grep("postpred", colnames(temp), value=TRUE)
      ]

      # * * * Compute p-values
      ICDM.p.values <- colMeans(temp.realized.values <= temp.postpred.values)


    } # closes if doing model checking


    # SUMMARIZING RESULTS -----------


    # Results from model fitting
    # Rename columns of results -------

    draws.from.MUPPET.model.blavaan.syntax.names <- draws.from.MUPPET.model

    # Draws from MUPPET model

    # Loop over parameters in the fitted MUPPET model
    # which.col=1
    # which.col=which.col+1
    for(which.col in 1:ncol(draws.from.MUPPET.model.blavaan.syntax.names)){

      # Proceed if this parameter is one of those that is not fixed
      # Basically, is it not fixed in the JAGS code (note, separate from supply values as data)

      if(colnames(draws.from.MUPPET.model.blavaan.syntax.names)[which.col] %in% combined.partable.jags$parameter.name.jags){

        # Get row from combined parameter table with this parameter
        temp <- combined.partable.jags %>% filter(parameter.name.jags == colnames(draws.from.MUPPET.model.blavaan.syntax.names)[which.col])

        # Paste together the blavaan syntax name as the new column name for the draws from the MUPPET model
        colnames(draws.from.MUPPET.model.blavaan.syntax.names)[which.col] <- paste0(temp$lhs, temp$op, temp$rhs)
      }
    }


    # Compute the summary statistics for the draws from the model -----
    draws.to.analyze <- draws.from.MUPPET.model.blavaan.syntax.names
    summary.statistics.MUPPET <- MCMCsummary(
      draws.to.analyze,
      HPD=TRUE,
      n.eff=FALSE,
      Rhat=FALSE,
      round=8,
      func=median,
      func_name = "median"
    )

    # Write out the summary statistics for the MUPPET model ------
    if(save.summary.stats.from.MUPPET){
      file.name=paste0("Summary Statistics MUPPET Model.csv")
      write.csv(
        x=summary.statistics.MUPPET,
        file=file.name
      )
    } # closes if saving summary statistics from MUPPET model


    # Write out draws from the MUPPET model -------
    if(save.draws.from.MUPPET){
      file.name <- "Draws from MUPPET model.out"
      to.write <- draws.from.MUPPET.model.blavaan.syntax.names
      write.table(to.write, file.name, row.names=FALSE)
    } # closes if saving draws from MUPPET model

    # Standardized solution to measurement model -----
    if(obtain.standardized.measurement.model){
      # Need to run twice to replace all the names

      # for(j in J:1){
      #   colnames(standardized.posterior.measurement.model) <- str_replace(colnames(standardized.posterior.measurement.model), paste0("x",j), indicators.names[j])
      # } # close loop over indicators
      # for(j in J:1){
      #   colnames(standardized.posterior.measurement.model) <- str_replace(colnames(standardized.posterior.measurement.model), paste0("x",j), indicators.names[j])
      # } # close loop over indicators

      # * Compute the summary statistics for standardized solution to the measurement model -----
      draws.to.analyze <- standardized.posterior.measurement.model
      summary.statistics.standardized.measurement <- MCMCsummary(
        draws.to.analyze,
        HPD=TRUE,
        n.eff=FALSE,
        Rhat=FALSE,
        round=8,
        func=median,
        func_name = "median"
      )

      # * Write out the summary statistics for the MUPPET model ------
      if(save.summary.stats.from.MUPPET){
        file.name=paste0("Summary Statistics Standardized Measurement Model.csv")
        write.csv(
          x=summary.statistics.standardized.measurement,
          file=file.name
        )
      } # closes if saving summary statistics from MUPPET model


    } # closes if obtaining standardized solution for the measurement model



    # Standardized solution to combined model -----
    if(obtain.standardized.combined.model){

      # * Compute standardized solution for the combined model ------
      standardized.posterior.combined.model <- standardizedPosterior.replace.draws(
        fitted.model.bsem.jags, replacement.draws = draws.from.MUPPET.model
      )

      # * For models with endogenous latent variables, replace standardized solution for measurement model params with that from the measurement model ------
      # The issue is that we're not enforcing the value for the marginal variance in models with exogenous latent variables
      # So when standarizing the results for the combined model, the values are slightly different than in the measurement model
      # So take the results for the measurement model

      if(combined.model.has.endogenous.latent){
        # Get names of observed variables from the measurement model

        ov.names.in.measurement.model <- lavNames(fitted.model.bcfa, type="ov")

        which.var=1
        # go through each variable in the measurement model
        for(which.var in 1:(length(ov.names.in.measurement.model))){
          which.col=1

          # go through each column in the results for the standardized posterior for the combined model
          for(which.col in 1:ncol(standardized.posterior.combined.model)){

            # Check if that column of the results for the standardized posterior for the combined model involves that variable from the measurement model
            if(grepl(ov.names.in.measurement.model[which.var], colnames(standardized.posterior.combined.model)[which.col])){

              # if it does, replace it with the column of the same name from the standardized solution to the measurement model
              standardized.posterior.combined.model[, colnames(standardized.posterior.combined.model)[which.col]] <- standardized.posterior.measurement.model[, colnames(standardized.posterior.combined.model)[which.col]]

            }
          }
        }

      } # closes if combined.model.has.endogenous.latent

      # * Write out standardized solution draws from the MUPPET model -------
      if(save.draws.from.MUPPET){
        file.name <- "Draws from MUPPET model standardized solution.out"
        to.write <- standardized.posterior.combined.model
        write.table(to.write, file.name, row.names=FALSE)
      } # closes if saving draws from MUPPET model, here for the standardized solution

      # * Compute the summary statistics for standardized solution to the measurement model -----
      draws.to.analyze <- standardized.posterior.combined.model
      summary.statistics.standardized.combined.model <- MCMCsummary(
        draws.to.analyze,
        HPD=TRUE,
        n.eff=FALSE,
        Rhat=FALSE,
        round=8,
        func=median,
        func_name = "median"
      )

      # * Write out the summary statistics for the standardized solution for the MUPPET model ------
      if(save.summary.stats.from.MUPPET){
        file.name=paste0("Summary Statistics Standardized Combined Model.csv")
        write.csv(
          x=summary.statistics.standardized.combined.model,
          file=file.name
        )
      } # closes if saving summary statistics from MUPPET model


    } # closes if obtaining standardized solution for the combined model





    # Results from model checking
    if(model.check){

      # * Rename columns of results -------

      # Need to run twice to replace all
      # Indicators
      for(j in J:1){
        colnames(PPMC.results.from.parallel) <- str_replace(colnames(PPMC.results.from.parallel), paste0("x",j), indicators.names[j])
      } # close loop over indicators

      for(j in J:1){
        colnames(PPMC.results.from.parallel) <- str_replace(colnames(PPMC.results.from.parallel), paste0("x",j), indicators.names[j])
      } # close loop over indicators


      # covariates
      if(model.has.covariates){
        for(h in H:1){
          colnames(PPMC.results.from.parallel) <- str_replace(colnames(PPMC.results.from.parallel), paste0("w",h), covariates.names[h])
        } # close loop over covariates

        for(h in H:1){
          colnames(PPMC.results.from.parallel) <- str_replace(colnames(PPMC.results.from.parallel), paste0("w",h), covariates.names[h])
        } # close loop over covariates
      } # closes if model has covariates

      # outcomes
      if(model.has.outcomes){
        for(k in K:1){
          colnames(PPMC.results.from.parallel) <- str_replace(colnames(PPMC.results.from.parallel), paste0("y",k), outcomes.names[k])
        } # close loop over outcomes
        for(k in K:1){
          colnames(PPMC.results.from.parallel) <- str_replace(colnames(PPMC.results.from.parallel), paste0("y",k), outcomes.names[k])
        } # close loop over outcomes
      } # closes if model has outcomes

      # colnames(PPMC.results.from.parallel)

      # * Write out the results for the realized and posterior predicted values from PPMC as a text file ------

      if(save.draws.from.MUPPET){
        data.file.name <- paste0("PPMC.Results.dat")
        write.csv(
          x=PPMC.results.from.parallel,
          file=paste0(PPMC.folder, data.file.name),
          row.names=FALSE
          #col.names=TRUE
        )
      }

      # * Compute summary statistics from PPMC  ------
      summary.statistics.PPMC <- MCMCsummary(
        PPMC.results.from.parallel,
        #params = parameters.to.summarize,
        HPD=TRUE,
        Rhat=FALSE, # won't compute because had to combine to one chain
        n.eff=FALSE, # won't compute because had to combine to one chain
        round=8,
        func=median,
        func_name = "median"
      )

      # * Append the p-values  ------
      summary.statistics.PPMC <- dplyr::mutate(summary.statistics.PPMC, p.post=NA)
      # summary.statistics.PPMC <- dplyr::mutate(summary.statistics.PPMC, rowname=rownames(summary.statistics.PPMC))

      # LR
      summary.statistics.PPMC[grep("realized.LR", rownames(summary.statistics.PPMC)), ]$p.post <- LR.p.value

      # Residual Correlation
      summary.statistics.PPMC[grep("realized_resid_cor", rownames(summary.statistics.PPMC)), ]$p.post <- resid.cor.p.values
      # names(resid.cor.p.values)==rownames(summary.statistics.PPMC[grep("realized_resid_cor", rownames(summary.statistics.PPMC)), ])

      # ICDM
      summary.statistics.PPMC[grep("realized_ICDM", rownames(summary.statistics.PPMC)), ]$p.post <- ICDM.p.values
      # names(ICDM.p.values)==rownames(summary.statistics.PPMC[grep("realized_ICDM", rownames(summary.statistics.PPMC)), ])


      # * Write out the summary statistics ------
      if(save.summary.stats.from.MUPPET){

        file.name=paste0("Summary Statistics PPMC.csv")
        write.csv(
          x=summary.statistics.PPMC,
          file=paste0(PPMC.folder, file.name)
        )

        # Write out p-values ------
        file.name=paste0("LR PPP-values.out")
        write.table(
          x=LR.p.value,
          file=paste0(PPMC.folder, file.name),
          row.names = FALSE,
          col.names = FALSE
        )

        file.name=paste0("Resid Cor PPP-values.out")
        write.table(
          x=resid.cor.p.values,
          file=paste0(PPMC.folder, file.name),
          row.names = FALSE,
          col.names = TRUE
        )

        file.name=paste0("ICDM PPP-values.out")
        write.table(
          x=ICDM.p.values,
          file=paste0(PPMC.folder, file.name),
          row.names = FALSE,
          col.names = TRUE
        )
      }



    } # closes if doing model checking

    # ASSEMBLE RESULTS TO BE RETURNED BY THE FUNCTION ---------
    # Foundational results
    MUPPET.CFA.function.result <- list(
      summary.statistics.MUPPET=summary.statistics.MUPPET,
      draws.from.MUPPET.model=draws.from.MUPPET.model
    )
    # If standardized solution for measurement model is requested
    if(obtain.standardized.measurement.model){
      standardized.measurement.list <- list(
        summary.statistics.standardized.measurement=summary.statistics.standardized.measurement,
        standardized.posterior.measurement.model=standardized.posterior.measurement.model
      )
      MUPPET.CFA.function.result <- append(MUPPET.CFA.function.result, standardized.measurement.list)
    }

    # If standardized solution for MUPPET model is requested
    if(obtain.standardized.combined.model){
      standardized.combined.list <- list(
        summary.statistics.standardized.combined.model=summary.statistics.standardized.combined.model,
        standardized.draws.from.MUPPET.model <- standardized.posterior.combined.model
      )
      MUPPET.CFA.function.result <- append(MUPPET.CFA.function.result, standardized.combined.list)
    }

    # If model check is requested
    if(model.check){
      model.check.list <- list(
        summary.statistics.PPMC=summary.statistics.PPMC,
        PPMC.results.from.parallel=PPMC.results.from.parallel,
        LR.p.value=LR.p.value,
        resid.cor.p.values=resid.cor.p.values,
        ICDM.p.values=ICDM.p.values
      )
      MUPPET.CFA.function.result <- append(MUPPET.CFA.function.result, model.check.list)
    }

    # TO BE DELETED: Append blavaan object for purposes of ripping things out for standardized solution
    MUPPET.CFA.function.result <- append(MUPPET.CFA.function.result, fitted.model.bsem.jags)
    names(MUPPET.CFA.function.result)[length(MUPPET.CFA.function.result)] <- "bsem.obj.cm"

    # Return the result of the function
    return(MUPPET.CFA.function.result)

  } # closes MUPPET.CFA.function

