##' A plotting method for the R mcmc model
##'
##' @title plot_r_mcmc
##' @param mcmc_ob The object returned from running the mcmc
##' @param burnin the burnin length
##' @return NULL
##' @author Jukka Ranta
##' @export
##' @importFrom graphics par
##' @importFrom graphics plot
##' @importFrom grDevices dev.new
plot_r_mcmc <- function(mcmc_ob, burnin){
    dev.new()
    par(mfrow = c(mcmc_ob$var_a$ns, 1), mar = c(2.5, 5.5, 2.5, 2.5))
    for(i in 1 : mcmc_ob$var_a$ns){
        plot(mcmc_ob$var_a$phi[burnin:mcmc_ob$var_a$MCMC, i],
             xlab="",
             ylab=mcmc_ob$var_b$data$sourcenames[i],
             pch=16,
             cex=0.5,
             cex.lab=1.7
             )
    }
}
