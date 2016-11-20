##' BUGS-model with multinomial(q[Z],1)-likelihoods 
##' for each human isolate with unknown group indicator (Z) 
##'
##' @title initialize_bugs_ST
##' @param ob The object for data formatting
##' @return A list
##' @export
##' @author Jukka Ranta
initialize_bugs_ST <- function(ob){
    nt <- length(ob$STu)
    sources <- ob$sourcesST
    Nisolates <- sum(ob$humansST)
    beta <- rep(1/ob$ns, ob$ns) # Dir prior parameters for sources
    data <- list(IST = ob$IST, sourcesST = ob$sourcesST, ns = ob$ns, nt = nt,
                 beta = beta, Nisolates = Nisolates)
    parameters <- c("phi", "Z", "qST")
    inits <- function(){
        list(g0 = rgamma(data$ns, 1, 1),
             q0 = structure(.Data = rep(1, data$ns * data$nt),
                            .Dim = c(data$ns, data$nt)),
             Z = round(runif(data$Nisolates,
                             0.5 / data$ns,
                             1 + 0.5 / data$ns) * data$ns)
             )
    }
    other <- list(ns = ob$ns,
                  sourcenames = ob$sourcenames)
    result <- list(data = data,
                   parameters = parameters,
                   inits = inits,
                   other = other)
    return(result)
}
