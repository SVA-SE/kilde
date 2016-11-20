##' plot_population_attribution.kilde_bugsmcmc_ST
##'
##' @title plot_population_attribution.kilde_bugsmcmc_ST
##' @param x An object of class kilde_bugsmcmc_ST
##' @param burnin The burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_population_attribution.kilde_bugsmcmc_ST <- function(x, burnin){
    ns <- x$other$ns
    STu <- x$var_b$STu
    sourcesST <- x$var_b$sourcesST
    MCMC <- x$bugs_result$n.iter
    qST <- x$bugs_result$sims.list$qST
    phi <- x$bugs_result$sims.list$phi
    sourcenames <- x$var_b$sourcenames
    MCMC <- MCMC - burnin
    burnin <- 0
    plot_population_attribution_ST_internal(ns,
                                            STu,
                                            sourcesST,
                                            MCMC,
                                            qST,
                                            phi,
                                            burnin,
                                            sourcenames)
}
#' plot_population_attribution.kilde_rmcmc_ST
##'
##' @title plot_population_attribution.kilde_rmcmc_ST
##' @param x An object of class kilde_rmcmc_ST
##' @param burnin The burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_population_attribution.kilde_rmcmc_ST <- function(x, burnin){
    ns <- x$var_a$ns
    STu <- x$var_a$STu
    sourcesST <- x$var_b$sourcesST
    MCMC <- x$var_a$MCMC
    qST <- x$var_a$qST
    phi <- x$var_a$phi
    sourcenames <- x$var_b$sourcenames
    plot_population_attribution_ST_internal(ns,
                                            STu,
                                            sourcesST,
                                            MCMC,
                                            qST,
                                            phi,
                                            burnin,
                                            sourcenames)
}
##' method tor plotting the ST model from R plot_population_attribution_ST_internal
##'
##' @title plot_population_attribution_ST_internal
##' @param ns ns 
##' @param STu STu 
##' @param sourcesST sourcesST 
##' @param MCMC MCMC 
##' @param phi phi 
##' @param burnin burnin 
##' @param sourcenames sourcenames 
##' @return A plot
##' @author Thomas Rosendal
plot_population_attribution_ST_internal <- function(ns,
                                                    STu,
                                                    sourcesST,
                                                    MCMC,
                                                    qST,
                                                    phi,
                                                    burnin,
                                                    sourcenames
                                                    ){
    ## The population attribution plot and the fit
    loST<-matrix(0, ns, length(STu))
    upST<-matrix(0, ns, length(STu))
                                        # check the fit of estimated freqs and observed freqs:
                                        # Set the number of picture frames as suitable:
    savepar <- par(mfrow=c(2, 1))
    on.exit(par(savepar))
    errorST <- matrix(NA, ns-1, length(STu))
    errorSThum <- numeric()
    mqST <- matrix(NA, ns, length(STu))
    fST <- numeric()
    for(i in 1:(ns - 1)){
        for(j in 1:length(STu)){
            fST[j] <- sourcesST[i, j] / sum(sourcesST[i, ])
        }
    }
    for(i in 1:ns){
        for(j in 1:length(STu)){
            mqST[i,j] <- mean(qST[burnin:MCMC, i, j]) 
            loup <- quantile(qST[burnin:MCMC, i, j],
                             c(0.025, 0.975),
                             names=FALSE)
            loST[i,j] <- loup[1]
            upST[i,j] <- loup[2]
        }
    }
    limx <- max(upST) * 1.01
    limy <- max(fST) * 1.1 
    plot(c(0, 1),
         c(0, 1),
         type = "l",
         xlim = c(0, limx),
         ylim = c(0, limy),
         xlab = "Estimated freq",
         ylab = "Observed sample freq",
         main = "ST")
    for(i in 1:(ns - 1)){
        errorST[i, 1:length(STu)] <- ( mqST[i, ] - sourcesST[i, ] / sum(sourcesST[i, ]))^2
        points(mqST[i, ],
               sourcesST[i, ] / sum(sourcesST[i, ]),
               col=i,
               pch=16)
        for(j in 1:length(STu)){
            points(c(loST[i, j], upST[i, j]),
                   rep(sourcesST[i, j] / sum(sourcesST[i, ]), 2),
                   col = i,
                   type = "l")}
    }
    ## start plotting SA results:
    s <- numeric()
    for(i in 1:ns){
        s[i] <- sd(phi[burnin:MCMC, i])
    }
    iix <- sort(s, index.return=TRUE)
    attach(iix)
    plot(density(phi[burnin:MCMC, ix[1]]),
         ylim = c(0, 60),
         col = ix[1],
         xlim = c(0, 1),
         xlab = "Source attribution",
         ylab = "", main = "", lwd = 2)
    for(i in 2:ns){
        points(density(phi[burnin:MCMC, ix[i]]),
               type = "l",
               col = ix[i],
               lwd = 2)
    }
    names <- c(sourcenames, rep("Unknown",
                                10))
    legend(x="topright",
           names[ix],
           text.col = ix)
}

errorsum <- numeric() 
for(i in 1:ns-1){
errorsum[i] <- sum(errorST[i,])
}
mphi <- numeric()
meanSThum <- numeric()

for(i in 1:ns){mphi[i] <- mean(phi[burnin:MCMC,i]) } 
for(j in 1:length(STu)){ meanSThum[j] <- sum(mphi*mqST[,j]);errorSThum[j]<-(humansST[j]/sum(humansST)-meanSThum[j])^2 }

ES <- sum(errorsum)   # sum of all squared errors
EShum <- sum(errorSThum)

