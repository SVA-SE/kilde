##' A method to run a model in bugs
##' @title run_model 
##' @param model The name of the model
##' This is the file name in the inst/models directory of the package
##' @param data A dataset
##' @param inits An inits function
##' @param parameters Parameters to monitor
##' @param n.chains The number of chains
##' @param n.burnin The burning length
##' @param n.iter The number of iterations
##' @return A model object
##' @author Thomas Rosendal
##' @export
run_model <- function(model = "SA_allelemodel1.jag",
                      df,
                      inits,
                      parameters,
                      n.chains = 1,
                      n.burnin = 99,
                      n.iter = 999) {
    bugs(df$bugs_data,
         inits,
         parameters,
         system.file(package = "kilde", paste0("models/", model)),
         n.chains = n.chains,
         n.burnin = n.burnin,
         n.iter = n.iter)
}
##' A method to run a model in jags
##' @title run_model_jags 
##' @param model The name of the model
##' This is the file name in the inst/models directory of the package
##' @param data A dataset
##' @param inits An inits function
##' @param parameters Parameters to monitor
##' @param n.chains The number of chains
##' @param n.burnin The burning length
##' @param n.iter The number of iterations
##' @return A model object
##' @author Thomas Rosendal
##' @export
run_model_jags<- function(model = "SA_allelemodel1.jag",
                      data,
                      inits,
                      parameters,
                      n.chains = 1,
                      n.burnin = 99,
                      n.iter = 999) {
    jags <- jags.model(system.file(package = "kilde",
                                   paste0("models/", model)),
                       data,
                       inits,
                       n.chains = n.chains,
                       n.adapt = n.burnin)                       
    update(jags, n.iter)
    jags.samples(jags, parameters, 1000)
}
