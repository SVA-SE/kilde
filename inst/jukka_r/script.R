DATA <- read.table("NMDD2015data_vers21_modif.txt",header=TRUE)
 
##attach(DATA)
z <- !is.na(DATA$ST)
 
Ctr = "FI"    # choose the Country  (DK, SE, NO, FI)
MCMC = 200   # choose number of MCMC iterations
burnin = 100  # choose burnin period
FULL = 0      # choose 1 for full Bayesian model (semi-supervised), but
              # choose 0 for separated estimation (supervised) of relative type frequencies
UM = 2        # choose 0 for setting 0 as 'data sample' for unknown source. 
              # choose 1 for taking the mean type frequencies as 'data sample' for unknown source.
              # choose 2 for taking the type frequencies drawn from the STs that were unique to humans, as 'data sample' for unknown source.             
source("dataformatting.R")
 
## ·         After this formatting is completed, the data should be
## ·         ready for either running the Gibbs sampling solely in R,
## ·         or running the corresponding BUGS model.  The MCMC
## ·         iterations (2100) works at least for FI data, but for SE
## ·         data it may be too high. These models are saving the MCMC
## ·         output for all parameters, which can take a lot of memory
## ·         if we have lots of isolates and lots of types and many
## ·         source populations. If it crashes, try with smaller MCMC
## ·         iterations.
 
## ·         After formatting, the MCMC running solely in R is done with the commands:
 
source("initializemcmc.R")
source("runmcmc.R")
source("plotting.R")
 
## ·         And the MCMC with the BUGS model:
 
library("R2OpenBUGS")
 
source("initializebugs.R")
 
# Below, BUGS model cannot handle large number of MCMC iterations for all parameters.
# Therefore, advisable to try with smaller number:
 
MCMC <- 200
burnin <- 100
 
res <- bugs(data,inits,parameters,"SA_allele_isolate.txt",n.chains=1,n.burnin=burnin,n.iter=MCMC)
 
attach.bugs(res)
 
# BUGS has returned iterations from burnin+1
 
MCMC <- MCMC-burnin
burnin <- 1
source("plotting.R")
