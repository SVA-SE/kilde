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
                 
                 runif(data$Nisolates, 0.5 / data$ns, 1 + 0.5 / data$ns) * data$ns),
             qASP = structure(.Data = rep(1 / data$nat[1], data$ns *
                                                           data$nat[1]),
                              .Dim = c(data$ns, data$nat[1])),
             qGLN = structure(.Data = rep(1 / data$nat[2], data$ns * data$nat[2]), .Dim = c(data$ns, data$nat[2])),
             qGLT = structure(.Data = rep(1 / data$nat[3], data$ns * data$nat[3]), .Dim = c(data$ns, data$nat[3])),
             qGLY = structure(.Data = rep(1 / data$nat[4], data$ns * data$nat[4]), .Dim = c(data$ns, data$nat[4])),
             qPGM = structure(.Data = rep(1 / data$nat[5], data$ns * data$nat[5]), .Dim = c(data$ns, data$nat[5])),
             qTKT = structure(.Data = rep(1 / data$nat[6], data$ns * data$nat[6]), .Dim = c(data$ns, data$nat[6])),
             qUNC = structure(.Data = rep(1 / data$nat[7], data$ns * data$nat[7]), .Dim = c(data$ns, data$nat[7])))
    }
    result <- list(data = data,
                   parameters = parameters,
                   inits = inits)
}
