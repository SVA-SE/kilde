
# you should have the original data table in your current working directory
# from which you read it as follows:

DATA <- read.table("NMDD2015data_vers21_modif.txt",header=TRUE)

attach(DATA)
z <- !is.na(ST)

Ctr = "FI"    # choose the Country  (DK, SE, NO, FI)
MCMC = 2100   # choose number of MCMC iterations
burnin = 100  # choose burnin period
FULL = 0      # choose 1 for full Bayesian model (semi-supervised), but 
              # choose 0 for separated estimation (supervised) of relative type frequencies 
UM = 2        # choose 0 for setting 0 as 'data sample' for unknown source.  
              # choose 1 for taking the mean type frequencies as 'data sample' for unknown source.
              # choose 2 for taking the type frequencies drawn from the STs that were unique to humans, as 'data sample' for unknown source.

              

source("dataformatting.R")

source("initializemcmc.R")

source("runmcmc.R")

source("plotting.R")

###############################
#BUGS model:

library("R2OpenBUGS")

source("initializebugs.R")

# Below, BUGS model cannot handle large number of MCMC iterations for all parameters.
# Therefore, advisable to try with smaller number:

MCMC <- 2000
burnin <- 100

res <- bugs(data,inits,parameters,"SA_allele_isolate.txt",n.chains=1,n.burnin=burnin,n.iter=MCMC)

attach.bugs(res)

# BUGS has returned iterations from burnin+1

MCMC <- MCMC-burnin
burnin <- 1
source("plotting.R")









