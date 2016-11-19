##' A plotting method for showing the mcmc history of the model
##'
##' @title plot_history
##' @param ns the number of sources
##' @param phi The phi object from the mcmc
##' @param MCMC The number of iterations
##' @param sourcenames A character vector of the sourcenames
##' @param burnin the burnin length
##' @return NULL
##' @author Jukka Ranta
##' @importFrom graphics par
##' @importFrom graphics plot
##' @importFrom grDevices dev.new
plot_history_internal <- function(ns,
                                  phi,
                                  MCMC,
                                  sourcenames,
                                  burnin) {
    parsave <- par(mfrow = c(ns, 1), mar = c(2.5, 5.5, 2.5, 2.5))
    on.exit(par(parsave))
    for(i in 1 : ns){
        plot(phi[burnin:MCMC, i],
             xlab = "",
             ylab = sourcenames[i],
             pch = 16,
             cex = 0.5,
             cex.lab = 1.7
             )
    }
}
##' plot_history
##' 
##' Defining the generic function
##' @title plot_history 
##' @param x The object to be plotted
##' @param burnin the burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_history <- function(x, burnin){
  UseMethod('plot_history')
}
plot_history.default = plot_history
##' History plotting function for kilde_rmcmc object
##'
##' 
##' @title plot_history.kilde_rmcmc
##' @param x The kilde_rmcmc class object
##' @param burnin The burning length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_history.kilde_rmcmc <- function(x, burnin){
    ns <- x$var_a$ns
    phi <- x$var_a$phi
    MCMC <- x$var_a$MCMC
    sourcenames <- x$var_b$data$sourcenames
    plot_history_internal(ns, phi, MCMC, sourcenames, burnin)
}
##' plot_history.kilde_bugsmcmc
##'
##' A plotting method for the history of the bugs object
##' 
##' @title plot_history.kilde_bugsmcmc
##' @param x The kilde_bugsmcmc class object
##' @param burnin the burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_history.kilde_bugsmcmc <- function(x, burnin){
    ns <- x$other$ns
    phi <- x$bugs_result$sims.list$phi
    MCMC <- x$bugs_result$n.iter - burnin
    sourcenames <- x$other$sourcenames
    burnin <- 0
    plot_history_internal(ns, phi, MCMC, sourcenames, burnin)
}
##' plot_population_attribution_internal: A plot of population attribution
##'
##' @title plot_population_attribution_internal
##' @param ns The number of sources
##' @param nat nat
##' @param MCMC The number of iterations 
##' @param phi The population attribution 
##' @param sourcenames A character vector of source names
##' @param burnin the burnin length
##' @param qASP qASP 
##' @param qGLN qGLN 
##' @param qGLT qGLT 
##' @param qGLY qGLY 
##' @param qPGM qPGM 
##' @param qTKT qTKT 
##' @param qUNC qUNC 
##' @param sourcesASP sourcesASP 
##' @param sourcesGLN sourcesGLN 
##' @param sourcesGLT sourcesGLT 
##' @param sourcesGLY sourcesGLY 
##' @param sourcesPGM sourcesPGM 
##' @param sourcesTKT sourcesTKT 
##' @param sourcesUNC sourcesUNC 
##' @return A plot
##' @importFrom graphics legend
##' @importFrom graphics points
##' @importFrom stats density
##' @importFrom stats quantile
##' @importFrom stats sd
##' @author Jukka Ranta
plot_population_attribution_internal <- function(ns,
                                                 nat,
                                                 MCMC,
                                                 phi,
                                                 sourcenames,
                                                 burnin,
                                                 qASP,
                                                 qGLN,
                                                 qGLT,
                                                 qGLY,
                                                 qPGM,
                                                 qTKT,
                                                 qUNC,
                                                 sourcesASP,
                                                 sourcesGLN,
                                                 sourcesGLT,
                                                 sourcesGLY,
                                                 sourcesPGM,
                                                 sourcesTKT,
                                                 sourcesUNC) {
    loASP <- matrix(0, ns, nat[1])
    loGLN <- matrix(0, ns, nat[2])
    loGLT <- matrix(0, ns, nat[3])
    loGLY <- matrix(0, ns, nat[4])
    loPGM <- matrix(0, ns, nat[5])
    loTKT <- matrix(0, ns, nat[6])
    loUNC <- matrix(0, ns, nat[7])
    ##
    upASP <- matrix(0, ns, nat[1])
    upGLN <- matrix(0, ns, nat[2])
    upGLT <- matrix(0, ns, nat[3])
    upGLY <- matrix(0, ns, nat[4])
    upPGM <- matrix(0, ns, nat[5])
    upTKT <- matrix(0, ns, nat[6])
    upUNC <- matrix(0, ns, nat[7])
    ##
    ## check the fit of estimated freqs and observed freqs:
    ## Set the number of picture frames as suitable:
    parsave <- par(mfrow=c(2,4))
    on.exit(par(parsave))
    errorASP <- matrix(NA, ns - 1, nat[1])
    errorGLN <- matrix(NA, ns - 1, nat[2])
    errorGLT <- matrix(NA, ns - 1, nat[3])
    errorGLY <- matrix(NA, ns - 1, nat[4])
    errorPGM <- matrix(NA, ns - 1, nat[5])
    errorTKT <- matrix(NA, ns - 1, nat[6])
    errorUNC <- matrix(NA, ns - 1, nat[7])
    errorASPhum <- numeric()
    errorGLNhum <- numeric()
    errorGLThum <- numeric()
    errorGLYhum <- numeric()
    errorPGMhum <- numeric()
    errorTKThum <- numeric()
    errorUNChum <- numeric()
    ## ASP plot
    mqASP <- matrix(NA, ns, nat[1])
    for(i in 1:ns){
        for(j in 1:nat[1]){
            mqASP[i, j] <- mean(qASP[burnin:MCMC, i, j]) 
            loup <- quantile(qASP[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loASP[i, j] <- loup[1]; upASP[i,j] <- loup[2]
        }
    }
    plot(c(0, 1),
         c(0, 1),
         type = "l",
         xlim = c(0, 1),
         ylim = c(0, 1),
         xlab = "Estimated freq",
         ylab = "Observed sample freq",
         main = "ASP")
    for(i in 1:(ns - 1)){
        errorASP[i, 1:nat[1]] <- (mqASP[i, ] - sourcesASP[i, ] /
                                                      sum(sourcesASP[i, ]))^2
        points(mqASP[i, ], sourcesASP[i, ] /
                           sum(sourcesASP[i, ]), col=i, pch=16)
        for(j in 1:nat[1]){ 
            points(c(loASP[i, j], upASP[i, j]),
                   rep(sourcesASP[i, j] /
                       sum(sourcesASP[i, ]), 2), col=i, type="l")
        }
    }
    ## GLN plot
    mqGLN <- matrix(NA, ns, nat[2])
    for(i in 1:ns){
        for(j in 1:nat[2]){
            mqGLN[i, j] <- mean(qGLN[burnin:MCMC, i, j]) 
            loup <- quantile(qGLN[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loGLN[i, j] <- loup[1]; upGLN[i,j] <- loup[2]
        }
    }
    plot(c(0, 1),
         c(0, 1),
         type = "l",
         xlim = c(0, 1),
         ylim = c(0, 1),
         xlab = "Estimated freq",
         ylab = "Observed sample freq",
         main = "GLN")
    for(i in 1:(ns - 1)){
        errorGLN[i, 1:nat[2]] <- (mqGLN[i, ] - sourcesGLN[i, ] /
                                                      sum(sourcesGLN[i, ]))^2
        points(mqGLN[i, ], sourcesGLN[i, ] /
                           sum(sourcesGLN[i, ]), col=i, pch=16)
        for(j in 1:nat[2]){ 
            points(c(loGLN[i, j], upGLN[i, j]),
                   rep(sourcesGLN[i, j] /
                       sum(sourcesGLN[i, ]), 2), col=i, type="l")
        }
    }
    ## GLT plot
    mqGLT <- matrix(NA, ns, nat[3])
    for(i in 1:ns){
        for(j in 1:nat[3]){
            mqGLT[i, j] <- mean(qGLT[burnin:MCMC, i, j]) 
            loup <- quantile(qGLT[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loGLT[i, j] <- loup[1]; upGLT[i,j] <- loup[2]
        }
    }
    plot(c(0, 1),
         c(0, 1),
         type = "l",
         xlim = c(0, 1),
         ylim = c(0, 1),
         xlab = "Estimated freq",
         ylab = "Observed sample freq",
         main = "GLT")
    for(i in 1:(ns - 1)){
        errorGLT[i, 1:nat[3]] <- (mqGLT[i, ] - sourcesGLT[i, ] /
                                                      sum(sourcesGLT[i, ]))^2
        points(mqGLT[i, ], sourcesGLT[i, ] /
                           sum(sourcesGLT[i, ]), col=i, pch=16)
        for(j in 1:nat[3]){ 
            points(c(loGLT[i, j], upGLT[i, j]),
                   rep(sourcesGLT[i, j] /
                       sum(sourcesGLT[i, ]), 2), col=i, type="l")
        }
    }
    ## GLY plot
    mqGLY <- matrix(NA, ns, nat[4])
    for(i in 1:ns){
        for(j in 1:nat[4]){
            mqGLY[i, j] <- mean(qGLY[burnin:MCMC, i, j]) 
            loup <- quantile(qGLY[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loGLY[i, j] <- loup[1]; upGLY[i,j] <- loup[2]
        }
    }
    plot(c(0, 1),
         c(0, 1),
         type = "l",
         xlim = c(0, 1),
         ylim = c(0, 1),
         xlab = "Estimated freq",
         ylab = "Observed sample freq",
         main = "GLY")
    for(i in 1:(ns - 1)){
        errorGLY[i, 1:nat[4]] <- (mqGLY[i, ] - sourcesGLY[i, ] /
                                                      sum(sourcesGLY[i, ]))^2
        points(mqGLY[i, ], sourcesGLY[i, ] /
                           sum(sourcesGLY[i, ]), col=i, pch=16)
        for(j in 1:nat[4]){ 
            points(c(loGLY[i, j], upGLY[i, j]),
                   rep(sourcesGLY[i, j] /
                       sum(sourcesGLY[i, ]), 2), col=i, type="l")
        }
    }
    ## PGM plot
    mqPGM <- matrix(NA, ns, nat[5])
    for(i in 1:ns){
        for(j in 1:nat[5]){
            mqPGM[i, j] <- mean(qPGM[burnin:MCMC, i, j]) 
            loup <- quantile(qPGM[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loPGM[i, j] <- loup[1]; upPGM[i,j] <- loup[2]
        }
    }
    plot(c(0, 1),
         c(0, 1),
         type = "l",
         xlim = c(0, 1),
         ylim = c(0, 1),
         xlab = "Estimated freq",
         ylab = "Observed sample freq",
         main = "PGM")
    for(i in 1:(ns - 1)){
        errorPGM[i, 1:nat[5]] <- (mqPGM[i, ] - sourcesPGM[i, ] /
                                                      sum(sourcesPGM[i, ]))^2
        points(mqPGM[i, ], sourcesPGM[i, ] /
                           sum(sourcesPGM[i, ]), col=i, pch=16)
        for(j in 1:nat[5]){ 
            points(c(loPGM[i, j], upPGM[i, j]),
                   rep(sourcesPGM[i, j] /
                       sum(sourcesPGM[i, ]), 2), col=i, type="l")
        }
    }
    ## TKT plot
    mqTKT <- matrix(NA, ns, nat[6])
    for(i in 1:ns){
        for(j in 1:nat[6]){
            mqTKT[i, j] <- mean(qTKT[burnin:MCMC, i, j]) 
            loup <- quantile(qTKT[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loTKT[i, j] <- loup[1]; upTKT[i,j] <- loup[2]
        }
    }
    plot(c(0, 1),
         c(0, 1),
         type = "l",
         xlim = c(0, 1),
         ylim = c(0, 1),
         xlab = "Estimated freq",
         ylab = "Observed sample freq",
         main = "TKT")
    for(i in 1:(ns - 1)){
        errorTKT[i, 1:nat[6]] <- (mqTKT[i, ] - sourcesTKT[i, ] /
                                                      sum(sourcesTKT[i, ]))^2
        points(mqTKT[i, ], sourcesTKT[i, ] /
                           sum(sourcesTKT[i, ]), col=i, pch=16)
        for(j in 1:nat[6]){ 
            points(c(loTKT[i, j], upTKT[i, j]),
                   rep(sourcesTKT[i, j] /
                       sum(sourcesTKT[i, ]), 2), col=i, type="l")
        }
    }
    ## UNC plot
    mqUNC <- matrix(NA, ns, nat[7])
    for(i in 1:ns){
        for(j in 1:nat[7]){
            mqUNC[i, j] <- mean(qUNC[burnin:MCMC, i, j]) 
            loup <- quantile(qUNC[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loUNC[i, j] <- loup[1]; upUNC[i,j] <- loup[2]
        }
    }
    plot(c(0, 1),
         c(0, 1),
         type = "l",
         xlim = c(0, 1),
         ylim = c(0, 1),
         xlab = "Estimated freq",
         ylab = "Observed sample freq",
         main = "UNC")
    for(i in 1:(ns - 1)){
        errorUNC[i, 1:nat[7]] <- (mqUNC[i, ] - sourcesUNC[i, ] /
                                                      sum(sourcesUNC[i, ]))^2
        points(mqUNC[i, ], sourcesUNC[i, ] /
                           sum(sourcesUNC[i, ]), col=i, pch=16)
        for(j in 1:nat[7]){ 
            points(c(loUNC[i, j], upUNC[i, j]),
                   rep(sourcesUNC[i, j] /
                       sum(sourcesUNC[i, ]), 2), col=i, type="l")
        }
    }
    ## start plotting SA results:
    s <- numeric()
    for(i in 1:ns){
        s[i] <- sd(phi[burnin:MCMC, i])
    }
    iix <- sort(s, index.return = TRUE)
    ##pdf("SAresult.pdf")
    plot(density(phi[burnin:MCMC, iix$ix[1]]),
         col = iix$ix[1],
         xlim = c(0, 1),
         xlab = "Posterior",
         ylab = "",
         main = "",
         lwd = 2)
    for(i in 2:ns){
        points(density(phi[burnin:MCMC, iix$ix[i]]),
               type = "l",
               col = iix$ix[i],
               lwd = 2)
    }
    names <- c(sourcenames, rep("Unknown", 10))
    legend(x = "topright",
           names[iix$ix],
           text.col=iix$ix)
}
##' Method for plotting population attribution
##'
##' @title plot_population_attribution
##' @param x the mcmc object 
##' @param burnin the burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_population_attribution <- function(x, burnin){
  UseMethod('plot_population_attribution')
}
plot_history.default = plot_population_attribution
##' A population attribution plot for the kilde_rmcmc object
##'
##' @title plot_population_attribution.kilde_rmcmc
##' @param x the mcmc_object of class kilde_rmcmc
##' @param burnin the burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_population_attribution.kilde_rmcmc <- function(x,
                                                    burnin) {
    ns <- x$var_a$ns
    nat <- x$var_b$inits$nat
    MCMC <- x$var_a$MCMC
    phi <- x$var_a$phi
    sourcenames <- x$var_b$data$sourcenames
    qASP <- x$var_a$qASP
    qGLN <- x$var_a$qGLN
    qGLT <- x$var_a$qGLT
    qGLY <- x$var_a$qGLY
    qPGM <- x$var_a$qPGM
    qTKT <- x$var_a$qTKT
    qUNC <- x$var_a$qUNC
    sourcesASP <- x$var_b$data$sourcesASP
    sourcesGLN <- x$var_b$data$sourcesGLN
    sourcesGLT <- x$var_b$data$sourcesGLT
    sourcesGLY <- x$var_b$data$sourcesGLY
    sourcesPGM <- x$var_b$data$sourcesPGM
    sourcesTKT <- x$var_b$data$sourcesTKT
    sourcesUNC <- x$var_b$data$sourcesUNC
    plot_population_attribution_internal(ns,
                                         nat,
                                         MCMC,
                                         phi,
                                         sourcenames,
                                         burnin,
                                         qASP,
                                         qGLN,
                                         qGLT,
                                         qGLY,
                                         qPGM,
                                         qTKT,
                                         qUNC,
                                         sourcesASP,
                                         sourcesGLN,
                                         sourcesGLT,
                                         sourcesGLY,
                                         sourcesPGM,
                                         sourcesTKT,
                                         sourcesUNC)
}
