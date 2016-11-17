##' run_bugs()
##'
##' a function that passes the formatted data to bugs 
##' 
##' @title run_bugs 
##' @param result The result from the initialize_bugs() function
##' @param MCMC The number of iterations in the MCMC
##' @param n.burnin the burnin
##' @param FULL Choose 1 for full Bayesian model (semi-supervised),
##'     but Choose 0 for separated estimation (supervised) of relative
##'     type frequencies
##' @param model The path to the model file
##' @param n.chains the number of chains
##' @return A list
##' @export
##' @author Thomas Rosendal
run_bugs <- function(result,
                     MCMC,
                     n.burnin,
                     FULL,
                     model = "SA_allelemodel3.jag",
                     n.chains = 1) {
    if(file.exists(system.file(package = "kilde", paste0("models/", model)))){
        modelfile <- system.file(package = "kilde", paste0("models/", model))    
    }else if(file.exists(model)) {
        modelfile <- model
    }else {
        stop("The model argument must be either the name of one of the models in the package or the complete path to a modelfile elsewhere")
    }
    result$data <- c(result$data, FULL = FULL)
    res2 <- bugs(result$data,
                 result$inits,
                 result$parameters,
                 model.file = modelfile,
                 n.chains = n.chains,
                 n.burnin = n.burnin,
                 n.iter = MCMC)
}
