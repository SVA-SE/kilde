library(kilde)
## Read in a format data
df <- sample_data()
df <- df[df$country == "Canada",]
ob <- dataformatting(DATA = df,
                     UM = 2)
######################################
## Initialize and then run mcmc in R
######################################
result <- initialize_mcmc(ns = ob$inits$ns,
                          nat = ob$inits$nat,
                          MCMC = 100,
                          Nisolates = ob$inits$Nisolates)
mcmc_ob <- runmcmc(result, ob, MCMC = 100, h = 0, FULL = 0)
##  Plot the results of this model.
plot_history(mcmc_ob, 50)
plot_population_attribution(mcmc_ob, 50)
##
################################################
## Initialize and then run the model in OpenBugs
################################################
##
## Below, BUGS model cannot handle a large number of MCMC iterations
## for all parameters.  Therefore, advisable to try with smaller
## number of iterations to start, perhaps 1000
##
initial_result <- initialize_bugs(ob)
result_bugs <- run_bugs(result = initial_result,
                        ob = ob,
                        MCMC = 1000,
                        n.burnin = 100,
                        FULL = 0)
plot_history(result_bugs, 100)
plot_population_attribution(result_bugs, 100)
