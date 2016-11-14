library(kilde)

## In order to make the sample data work I needed to add some
## capitalization to a few strings
df <- sample_data()
df <- df[df$country == "Canada",]
ob <- dataformatting(DATA = df,
                     UM = 2)
result <- initialize_mcmc(ns = ob$inits$ns,
                          nat = ob$inits$nat,
                          MCMC = 100,
                          Nisolates = ob$inits$Nisolates)
mcmc_ob <- runmcmc(result, ob, MCMC = 100, h = 0, FULL = 0)
plot_r_mcmc(mcmc_ob, 50)

## ## ·         And the MCMC with the BUGS model:
 
## library("R2OpenBUGS")
 
## source("initializebugs.R")
 
## # Below, BUGS model cannot handle large number of MCMC iterations for all parameters.
## # Therefore, advisable to try with smaller number:
 
## MCMC <- 200
## burnin <- 100
 
## res <- bugs(data,inits,parameters,"SA_allele_isolate.txt",n.chains=1,n.burnin=burnin,n.iter=MCMC)
 
## attach.bugs(res)
 
## # BUGS has returned iterations from burnin+1
 
## MCMC <- MCMC-burnin
## burnin <- 1
## source("plotting.R")
