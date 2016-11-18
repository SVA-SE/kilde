
DATA <- read.table("NMDD2015data_vers21_modif.txt",header=TRUE)

attach(DATA)
# use ALL allele types
z <- !is.na(ST)

Ctr = "FI"     # choose the Country  (DK, SE, NO, FI)
MCMC <- 1000   # number of MCMC iterations
burnin <- 101  # set the burn-in period length
FULL <- 0      # choose 1 for full Bayesian model (semi-supervised), but 
               # choose 0 for separated estimation (supervised) of relative type frequencies 
UM = 2         # choose 0 for setting 0 as 'data sample' for unknown source.  
               # choose 1 for taking the mean type frequencies as 'data sample' for unknown source.
               # choose 2 for taking the type frequencies drawn from the STs that were unique to humans, as 'data sample' for unknown source.

source("dataformattingST.R")



source("initializemcmcST.R")
source("runmcmcST.R")
source("plottingST.R")


############################
#### try BUGS models:
library("R2OpenBUGS")

source("initializebugsST.R")

MCMC <- 1000
burnin <- 100

resST <- bugs(data,inits,parameters,"SA_ST_isolate.txt",n.chains=1,n.burnin=burnin,n.iter=MCMC)
attach.bugs(resST)

MCMC <- MCMC-burnin
burnin <- 1
source("plottingST.R")


