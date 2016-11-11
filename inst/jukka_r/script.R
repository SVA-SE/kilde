library(kilde)
DATA <- read.table("NMDD2015data_vers21_modif.txt",header=TRUE)
Ctr = "FI"    # choose the Country  (DK, SE, NO, FI)
UM = 2        # choose 0 for setting 0 as 'data sample' for unknown source.
              # choose 1 for taking the mean type frequencies as 'data sample' for unknown source.
              # choose 2 for taking the type frequencies drawn from the STs that were unique to humans, as 'data sample' for unknown source.             
##source("dataformatting.R")
ob <- dataformatting(DATA, Ctr, UM)

source("initializemcmc.R")
result <- initialize_mcmc(ob$inits$ns, ob$inits$nat, MCMC, ob$inits$Nisolates)

source("runmcmc.R")
mcmc_ob <- runmcmc(result, ob, MCMC = 150, h = 0, FULL = 0)

source("plotting.R")
plota(mcmc_ob, 100)

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
