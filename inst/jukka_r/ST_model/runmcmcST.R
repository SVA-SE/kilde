##' A method to run the mcmc in R: ST model
##'
##' @title runmcmc_ST
##' @param result An object from the initialize_mcmc_ST function
##' @param ob An object from the dataformatting function
##' @param h A control parameter (for avoiding very small values of q
##'     ~Dirich, maybe never needed)
##' @param FULL Choose 1 for full Bayesian model (semi-supervised),
##'     but Choose 0 for separated estimation (supervised) of relative
##'     type frequencies
##' @return A list of objects to be passed onto the plot function
##' @export
##' @author Jukka Ranta
runmcmc_ST<- function(result,
                      ob,
                      h = 0,
                      FULL = 0){
    for(mc in 2:result$MCMC){ 
        ## FULL CONDITIONAL DISTRIBUTIONS FOR 
        ## THE RELATIVE ST TYPE FREQUENCIES IN EACH SOURCE POPULATION:
        for(i in 1:result$ns){   ## ns = number of source populations
            for(j in 1:length(result$STu)){ 
                result$STpar0[i, j] <- rgamma(1,
                                              ob$sourcesST[i, j] +
                                              FULL * sum(result$S[, i] * ob$IST[, j] > 0) +
                                              1 / length(result$STu),
                                              1)
                result$qST[mc, i, j] <- (result$STpar0[i, j] + h) / (sum(result$STpar0[i, ]) +
                                                                     h * length(result$STu))
            }
        }  ## end of ns
        ## FULL CONDITIONAL DISTRIBUTION FOR EACH Z: 
        ## (Z = group indicators for each isolate) 
        for(i in 1:result$Nisolates){
            for(k in 1:result$ns){ 
                result$phii0[i, k] <- exp(log(result$phi[mc - 1, k]) +
                                          log(result$qST[mc, k, ob$ind[i]])
                                          )
            }
            for(k in 1:result$ns){
                result$phii[i, k] <- exp(log(result$phii0[i, k]) -
                                         log(sum(result$phii0[i, ]))
                                         ) # normalize
            }
            nn <- 1:result$ns
            result$Z[mc, i] <- min(nn[runif(1) < cumsum(result$phii[i, ])])     ## "rcat(1,phii[i,])"
            
            for(k in 1:result$ns){
                result$S[i, k] <- result$Z[mc, i] == k   ## membership of ith isolate to group k  
            }
        } ## end Nisolates
        
        ## FULL CONDITIONAL DISTRIBUTION FOR phi:
        ## (phi = proportions of source populations)
        for(k in 1:result$ns){
            result$phi0[k] <- rgamma(1,
                                     sum(result$S[, k]) + 1 / result$ns,
                                     1)
        }
        for(k in 1:result$ns){
            result$phi[mc, k] <- result$phi0[k] / sum(result$phi0[])
        }   
    }  ## end of MCMC
    result <- list(var_a = result, var_b = ob)
    class(result) <- "kilde_rmcmc"
    return(result)
}
