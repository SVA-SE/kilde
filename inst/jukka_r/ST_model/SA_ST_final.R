library(kilde)
DATA <- read.table("../NMDD2015data_vers21_modif.txt",header=TRUE)
DATA <- DATA[DATA$Country == "FI",]
## First the ST model in R
ob <- dataformatting_ST(DATA = DATA, UM = 2)
result <- initialize_mcmc_ST(ob$ns, ob$STu, MCMC = 100, ob$Nisolates)
mcmc_ob <- runmcmc_ST(result = result,
           ob = ob,
           h = 0,
           FULL = 0)
plot_history(mcmc_ob, 50)
## Then the same model in BUGS
initial_result <- initialize_bugs_ST(ob)
result_bugs <- kilde::run_bugs(result = initial_result,
                               ob = ob,
                               MCMC = 1000,
                               n.burnin = 100,
                               FULL = 0,
                               model = "SA_ST_model.jag",
                               n.chains = 1)
plot_history(result_bugs, 100)



