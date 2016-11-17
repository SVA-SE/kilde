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
plot_r_mcmc(mcmc_ob, 50)
plot_posterior_mcmc(mcmc_ob, 50)
## 
## And the MCMC with the BUGS model: Below, BUGS model cannot handle
## large number of MCMC iterations for all parameters.  Therefore,
## advisable to try with smaller number:
##
## Initialize and then run the model in bugs
result <- initialize_bugs(ob)
result <- run_bugs(result = result,
                   MCMC = 1000,
                   n.burnin = 100,
                   FULL = 0)
## attach.bugs(res)
 
## # BUGS has returned iterations from burnin+1
 
## MCMC <- MCMC-burnin
## burnin <- 1
## source("plotting.R")
