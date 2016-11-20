library(kilde)
DATA <- read.table("../NMDD2015data_vers21_modif.txt",header=TRUE)
DATA <- DATA[DATA$Country == "FI",]

MCMC <- 1000   ## number of MCMC iterations
burnin <- 101  ## set the burn-in period length
FULL <- 0      ## choose 1 for full Bayesian model (semi-supervised), but 

source("dataformattingST.R")
ob <- dataformatting_ST(DATA = DATA, UM = 2)

source("initializemcmcST.R")
result <- initialize_mcmc_ST(ob$ns, ob$STu, MCMC = 100, ob$Nisolates)

source("runmcmcST.R")
mcmc_ob <- runmcmc_ST(result = result,
           ob = ob,
           h = 0,
           FULL = 0)
##debugonce(kilde:::plot_history.kilde_rmcmc)
plot_history(mcmc_ob, 50)

source("plottingST.R")

############################
#### try BUGS models:
source("initializebugsST.R")
initial_result <- initialize_bugs_ST(ob)
result_bugs <- kilde::run_bugs(result = initial_result,
                               ob = ob,
                               MCMC = 500,
                               n.burnin = 100,
                               FULL = 0,
                               model = "SA_ST_model.jag",
                               n.chains = 1)
plot_history(result_bugs, 100)

resST <- bugs(data,inits,parameters,"SA_ST_isolate.txt",n.chains=1,n.burnin=burnin,n.iter=MCMC)
attach.bugs(resST)

MCMC <- MCMC-burnin
burnin <- 1
source("plottingST.R")


