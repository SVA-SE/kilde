library(kilde)
## Read in a format data
df <- sample_data()
df <- df[df$country == "Canada",]
ob <- dataformatting(DATA = df,
                     UM = 2)
## Initialize and then run mcmc in R
result <- initialize_mcmc(ns = ob$inits$ns,
                          nat = ob$inits$nat,
                          MCMC = 100,
                          Nisolates = ob$inits$Nisolates)
mcmc_ob <- runmcmc(result, ob, MCMC = 100, h = 0, FULL = 0)
##  Plot it
plot_history(mcmc_ob, 50)
plot_population_attribution(mcmc_ob, 50)
## 
## And the MCMC with the BUGS model: Below, BUGS model cannot handle
## large number of MCMC iterations for all parameters.  Therefore,
## advisable to try with smaller number:
##
## Initialize and then run the model in bugs
initial_result <- initialize_bugs(ob)

result_bugs <- run_bugs(result = initial_result,
                        MCMC = 1000,
                        n.burnin = 100,
                        FULL = 0)
plot_history(result_bugs, 100)

## How to plot this?
## result$

## str(df)

library(R2OpenBUGS)
attach.bugs(result_bugs)
str(result_bugs$sims.list)
str(result_bugs$n.chains)
str(result_bugs$n.iter)
str(result_bugs$n.burnin)
plot_r_mcmc()


## ## BUGS has returned iterations from burnin+1

## MCMC <- MCMC-burnin
## burnin <- 1
## source("plotting.R")
