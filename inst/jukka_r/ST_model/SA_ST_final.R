DATA <- read.table("../NMDD2015data_vers21_modif.txt",header=TRUE)
DATA <- DATA[DATA$Country == "FI",]

MCMC <- 1000   ## number of MCMC iterations
burnin <- 101  ## set the burn-in period length
FULL <- 0      ## choose 1 for full Bayesian model (semi-supervised), but 

source("dataformattingST.R")
ob <- dataformatting_ST(DATA = DATA, UM = 2)

source("initializemcmcST.R")
result <- initialize_mcmc_ST(ob$ns, ob$STu, MCMC = 1000, ob$Nisolates)

source("runmcmcST.R")
source("plottingST.R")


############################
#### try BUGS models:
library("R2OpenBUGS")
source("initializebugsST.R")
result_bugs <- initialize_bugs_ST(ob)

MCMC <- 1000
burnin <- 100

resST <- bugs(data,inits,parameters,"SA_ST_isolate.txt",n.chains=1,n.burnin=burnin,n.iter=MCMC)
attach.bugs(resST)

MCMC <- MCMC-burnin
burnin <- 1
source("plottingST.R")


