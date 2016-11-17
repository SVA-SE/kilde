## VARIABLE STRUCTURES FOR THE MCMC (GIBBS) SAMPLER:
## This script requires a few variables in the global env
##' A method to initialize the data objects for mcmc
##'
##' @title initialize_mcmc
##' @param ns an integer, the number of sources
##' @param nat A vector the the number of alleles at each loci. In the
##'     Order: c("ASP", "GLN", "GLT", "GLY", "PGM", "TKT", "UNC")
##' @param MCMC the Number of iterations for the MCMC
##' @param Nisolates The number of human isolates
##' @return a list of the initial values for the simulation
##' @importFrom stats rgamma
##' @importFrom stats runif
##' @export
##' @author Jukka Ranta
initialize_mcmc <- function(ns, nat, MCMC, Nisolates){
    ASPpar0 <- matrix(0, ns, nat[1])
    GLNpar0 <- matrix(0, ns, nat[2])
    GLTpar0 <- matrix(0, ns, nat[3])
    GLYpar0 <- matrix(0, ns, nat[4])
    PGMpar0 <- matrix(0, ns, nat[5])
    TKTpar0 <- matrix(0, ns, nat[6])
    UNCpar0 <- matrix(0, ns, nat[7])
    qASP <- structure(.Data=c(rep(0, MCMC * ns * nat[1])), .Dim=c(MCMC, ns, nat[1]))
    qGLN <- structure(.Data=c(rep(0, MCMC * ns * nat[2])), .Dim=c(MCMC, ns, nat[2]))
    qGLT <- structure(.Data=c(rep(0, MCMC * ns * nat[3])), .Dim=c(MCMC, ns, nat[3]))
    qGLY <- structure(.Data=c(rep(0, MCMC * ns * nat[4])), .Dim=c(MCMC, ns, nat[4]))
    qPGM <- structure(.Data=c(rep(0, MCMC * ns * nat[5])), .Dim=c(MCMC, ns, nat[5]))
    qTKT <- structure(.Data=c(rep(0, MCMC * ns * nat[6])), .Dim=c(MCMC, ns, nat[6]))
    qUNC <- structure(.Data=c(rep(0, MCMC * ns * nat[7])), .Dim=c(MCMC, ns, nat[7]))
    phii <- matrix(0, Nisolates, ns)
    phii0 <- matrix(0, Nisolates, ns)
    Z <- matrix(0, MCMC, Nisolates)
    S <- matrix(0, Nisolates, ns) 
    phi0 <- numeric()
    phi <- matrix(0, MCMC, ns) 
    ## INITIAL VALUES FOR MCMC (GIBBS):
    for(i in 1 : ns){
        for(j in 1 : nat[1]){qASP[1, i, j] <- 1/nat[1]} 
        for(j in 1 : nat[2]){qGLN[1, i, j] <- 1/nat[2]} 
        for(j in 1 : nat[3]){qGLT[1, i, j] <- 1/nat[3]} 
        for(j in 1 : nat[4]){qGLY[1, i, j] <- 1/nat[4]} 
        for(j in 1 : nat[5]){qPGM[1, i, j] <- 1/nat[5]} 
        for(j in 1 : nat[6]){qTKT[1, i, j] <- 1/nat[6]} 
        for(j in 1 : nat[7]){qUNC[1, i, j] <- 1/nat[7]} 
    }
    ga <- rgamma(ns, 1, 1)
    phi[1, 1 : ns] <- ga/sum(ga)  ## source probabilities
    ## Indicators for group membership of human isolates:
    Z[1, 1 : Nisolates] <- round(runif(Nisolates, 0.5/ns, 1 + 0.5/ns) * ns)
    for(s in 1 : Nisolates){
        for(k in 1:ns){
            S[s, k] <- Z[1, s] == k   # membership of sth isolate to group k  
        }
    }
    result <- list(ASPpar0 = ASPpar0, GLNpar0 = GLNpar0, 
                   GLTpar0 = GLTpar0, GLYpar0 = GLYpar0, MCMC = MCMC, 
                   nat = nat, Nisolates = Nisolates, ns = ns, 
                   PGMpar0 = PGMpar0, phi = phi, phi0 = phi0, 
                   phii = phii, phii0 = phii0, qASP = qASP, 
                   qGLN = qGLN, qGLT = qGLT, qGLY = qGLY, qPGM = qPGM, 
                   qTKT = qTKT, qUNC = qUNC, S = S, TKTpar0 = TKTpar0, 
                   UNCpar0 = UNCpar0, Z = Z)
     return(result)   
}
##' initialize_bugs
##'
##' Packages a number of datasets from the data cleaning function to
##' prepare them for open bugs
##' @title initialize_bugs
##' @param ob the object from the data formatting function
##' @return An object that can be passed to the bug mcmc engine
##' @author Jukka Ranta
##' @export
initialize_bugs <- function(ob) {
    beta <- rep(1 / ob$inits$ns, ob$inits$ns)  # Dir prior parameters for sources
    data <- list(IASP = ob$data$IASP, IGLN = ob$data$IGLN,
                 IGLT = ob$data$IGLT, IGLY = ob$data$IGLY,
                 IPGM = ob$data$IPGM, ITKT = ob$data$ITKT,
                 IUNC = ob$data$IUNC, sourcesASP = ob$data$sourcesASP,
                 sourcesGLN = ob$data$sourcesGLN,
                 sourcesGLT = ob$data$sourcesGLT,
                 sourcesGLY = ob$data$sourcesGLY,
                 sourcesPGM = ob$data$sourcesPGM,
                 sourcesTKT = ob$data$sourcesTKT,
                 sourcesUNC = ob$data$sourcesUNC,  ns = ob$inits$ns,
                 nat = ob$inits$nat, beta = beta,
                 Nisolates = ob$inits$Nisolates)
    parameters <- c("qASP", "qGLN", "qGLT",
                    "qGLY", "qPGM", "qTKT",
                    "qUNC", "phi", "P", "Z")
    inits <- function(){
        list(g0 = rgamma(data$ns, 1, 1),
             Z = round(
                 runif(data$Nisolates, 0.5 / data$ns, 1 + 0.5 / data$ns) * data$ns
             ),
             qASP = structure(.Data = rep(1 / data$nat[1], data$ns * data$nat[1]),
                              .Dim = c(data$ns, data$nat[1])),
             qGLN = structure(.Data = rep(1 / data$nat[2], data$ns * data$nat[2]),
                              .Dim = c(data$ns, data$nat[2])),
             qGLT = structure(.Data = rep(1 / data$nat[3], data$ns * data$nat[3]),
                              .Dim = c(data$ns, data$nat[3])),
             qGLY = structure(.Data = rep(1 / data$nat[4], data$ns * data$nat[4]),
                              .Dim = c(data$ns, data$nat[4])),
             qPGM = structure(.Data = rep(1 / data$nat[5], data$ns * data$nat[5]),
                              .Dim = c(data$ns, data$nat[5])),
             qTKT = structure(.Data = rep(1 / data$nat[6], data$ns * data$nat[6]),
                              .Dim = c(data$ns, data$nat[6])),
             qUNC = structure(.Data = rep(1 / data$nat[7], data$ns * data$nat[7]),
                              .Dim = c(data$ns, data$nat[7])))
    }
    result <- list(data = data,
                   parameters = parameters,
                   inits = inits)
}
