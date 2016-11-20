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
##' History plotting function for kilde_rmcmc_ST object
##'
##' 
##' @title plot_history.kilde_rmcmc_ST
##' @param x The kilde_rmcmc_ST class object
##' @param burnin The burning length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_history.kilde_rmcmc_ST<- function(x, burnin){
    ns <- x$var_a$ns
    phi <- x$var_a$phi
    MCMC <- x$var_a$MCMC
    sourcenames <- x$var_b$sourcenames
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
##' plot_history.kilde_bugsmcmc_ST
##'
##' A plotting method for the history of the bugs object
##' 
##' @title plot_history.kilde_bugsmcmc_ST
##' @param x The kilde_bugsmcmc_ST class object
##' @param burnin the burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_history.kilde_bugsmcmc_ST<- function(x, burnin){
    ns <- x$other$ns
    phi <- x$bugs_result$sims.list$phi
    MCMC <- x$bugs_result$n.iter - burnin
    sourcenames <- x$other$sourcenames
    burnin <- 0
    plot_history_internal(ns, phi, MCMC, sourcenames, burnin)
}
##' plot_modelfit_internal: A plot of population attribution
##'
##' @title plot_modelfit_internal
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
plot_modelfit_internal <- function(ns,
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
##' @title plot_modelfit
##' @param x the mcmc object 
##' @param burnin the burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_modelfit <- function(x, burnin){
    UseMethod('plot_modelfit')
}
plot_modelfit.default = plot_modelfit
##' A population attribution plot for the kilde_rmcmc object
##'
##' @title plot_modelfit.kilde_rmcmc
##' @param x the mcmc_object of class kilde_rmcmc
##' @param burnin the burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_modelfit.kilde_rmcmc <- function(x,
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
    plot_modelfit_internal(ns,
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
##' plot_modelfit.kilde_bugsmcmc
##'
##' A method for ploting the population attribution from the bugs code
##'
##' @title plot_modelfit.kilde_bugsmcmc
##' @param x the kilde_bugsmcmc object
##' @param burnin The burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_modelfit.kilde_bugsmcmc <- function(x,
                                         burnin) {
    ## pick out objects to pass to plotting
    ns <- x$other$ns
    nat <- x$var_b$inits$nat
    MCMC <- x$bugs_result$n.iter
    phi <- x$bugs_result$sims.list$phi
    sourcenames <- x$other$sourcenames
    qASP <- x$bugs_result$sims.list$qASP
    qGLN <- x$bugs_result$sims.list$qGLN
    qGLT <- x$bugs_result$sims.list$qGLT
    qGLY <- x$bugs_result$sims.list$qGLY
    qPGM <- x$bugs_result$sims.list$qPGM
    qTKT <- x$bugs_result$sims.list$qTKT
    qUNC <- x$bugs_result$sims.list$qUNC
    sourcesASP <- x$var_b$data$sourcesASP
    sourcesGLN <- x$var_b$data$sourcesGLN
    sourcesGLT <- x$var_b$data$sourcesGLT
    sourcesGLY <- x$var_b$data$sourcesGLY
    sourcesPGM <- x$var_b$data$sourcesPGM
    sourcesTKT <- x$var_b$data$sourcesTKT
    sourcesUNC <- x$var_b$data$sourcesUNC
    MCMC <- MCMC - burnin
    burnin <- 0
    ## Run plotting function
    plot_modelfit_internal(ns,
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
##' The internal function to plot the sample attribution
##'
##' @title plot_sample_attribution_internal
##' @param Z Z
##' @param burnin the burning length
##' @param MCMC The number of iterations
##' @param ns the number of sources
##' @param sourcenames A character vector of the source names
##' @importFrom graphics hist
##' @return A plot
##' @author Jukka Ranta
plot_sample_attribution_internal <- function(Z,
                                             burnin,
                                             MCMC,
                                             ns,
                                             sourcenames){
    hist(Z[burnin:MCMC, ],
         0.5:(ns + 0.5),
         col = "turquoise3",
         labels = c(sourcenames, "Other"),
         axes = FALSE,
         ylab = "",
         xlab = "",
         main = "P(Z | data)")
}
##' Method for plotting sample attribution
##'
##' @title plot_sample_attribution
##' @return A plot
##' @export
##' @author Thomas Rosendal
##' @param x the model object
##' @param burnin the length of the burnin
plot_sample_attribution <- function(x, burnin){
    UseMethod('plot_sample_attribution')
}
plot_sample_attribution.default = plot_sample_attribution
##' plot_sample_attribution.kilde_rmcmc
##'
##' 
##' @title plot_sample_attribution.kilde_rmcmc
##' @param x the model output object
##' @param burnin the burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_sample_attribution.kilde_rmcmc <- function(x,
                                                burnin) {
    Z <- x$var_a$Z
    ns <- x$var_a$ns
    MCMC <- x$var_a$MCMC
    sourcenames <- x$var_b$data$sourcenames
    plot_sample_attribution_internal(Z,
                                     burnin,
                                     MCMC,
                                     ns,
                                     sourcenames)
}
##' plot_sample_attribution.kilde_rmcmc_ST
##'
##' 
##' @title plot_sample_attribution.kilde_rmcmc_ST
##' @param x the model output object
##' @param burnin the burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_sample_attribution.kilde_rmcmc_ST<- function(x,
                                                  burnin) {
    Z <- x$var_a$Z
    ns <- x$var_a$ns
    MCMC <- x$var_a$MCMC
    sourcenames <- x$var_b$sourcenames
    plot_sample_attribution_internal(Z,
                                     burnin,
                                     MCMC,
                                     ns,
                                     sourcenames)
}
##' plot_sample_attribution.kilde_bugsmcmc
##'
##' @title plot_sample_attribution.kilde_bugsmcmc
##' @param x the model object
##' @param burnin The burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_sample_attribution.kilde_bugsmcmc <- function(x,
                                                   burnin) {
    ## pick out objects to pass to plotting
    Z <- x$bugs_result$sims.list$Z
    ns <- x$other$ns
    MCMC <- x$bugs_result$n.iter
    sourcenames <- x$other$sourcenames
    MCMC <- MCMC - burnin
    burnin <- 0
    ## Run plotting function
    plot_sample_attribution_internal(Z,
                                     burnin,
                                     MCMC,
                                     ns,
                                     sourcenames)
}
##' plot_sample_attribution.kilde_bugsmcmc_ST
##'
##' @title plot_sample_attribution.kilde_bugsmcmc_ST
##' @param x the model object
##' @param burnin The burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_sample_attribution.kilde_bugsmcmc_ST<- function(x,
                                                     burnin) {
    ## pick out objects to pass to plotting
    Z <- x$bugs_result$sims.list$Z
    ns <- x$other$ns
    MCMC <- x$bugs_result$n.iter
    sourcenames <- x$other$sourcenames
    MCMC <- MCMC - burnin
    burnin <- 0
    ## Run plotting function
    plot_sample_attribution_internal(Z,
                                     burnin,
                                     MCMC,
                                     ns,
                                     sourcenames)
}
##' plot_modelfit.kilde_bugsmcmc_ST
##'
##' @title plot_modelfit.kilde_bugsmcmc_ST
##' @param x An object of class kilde_bugsmcmc_ST
##' @param burnin The burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_modelfit.kilde_bugsmcmc_ST <- function(x, burnin){
    ns <- x$other$ns
    STu <- x$var_b$STu
    sourcesST <- x$var_b$sourcesST
    MCMC <- x$bugs_result$n.iter
    qST <- x$bugs_result$sims.list$qST
    phi <- x$bugs_result$sims.list$phi
    sourcenames <- x$var_b$sourcenames
    MCMC <- MCMC - burnin
    burnin <- 0
    plot_modelfit_ST_internal(ns,
                              STu,
                              sourcesST,
                              MCMC,
                              qST,
                              phi,
                              burnin,
                              sourcenames)
}
##' plot_modelfit.kilde_rmcmc_ST
##'
##' @title plot_modelfit.kilde_rmcmc_ST
##' @param x An object of class kilde_rmcmc_ST
##' @param burnin The burnin length
##' @return A plot
##' @export
##' @author Thomas Rosendal
plot_modelfit.kilde_rmcmc_ST <- function(x, burnin){
    ns <- x$var_a$ns
    STu <- x$var_a$STu
    sourcesST <- x$var_b$sourcesST
    MCMC <- x$var_a$MCMC
    qST <- x$var_a$qST
    phi <- x$var_a$phi
    sourcenames <- x$var_b$sourcenames
    plot_modelfit_ST_internal(ns,
                              STu,
                              sourcesST,
                              MCMC,
                              qST,
                              phi,
                              burnin,
                              sourcenames)
}
##' method tor plotting the ST model from R plot_modelfit_ST_internal
##'
##' @title plot_modelfit_ST_internal
##' @param ns ns 
##' @param STu STu 
##' @param sourcesST sourcesST 
##' @param MCMC MCMC 
##' @param qST qST 
##' @param phi phi 
##' @param burnin burnin 
##' @param sourcenames sourcenames 
##' @return A plot
##' @author Thomas Rosendal
plot_modelfit_ST_internal <- function(ns,
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
    ##attach(iix)
    plot(density(phi[burnin:MCMC, iix$ix[1]]),
         ylim = c(0, 60),
         col = iix$ix[1],
         xlim = c(0, 1),
         xlab = "Source attribution",
         ylab = "", main = "", lwd = 2)
    for(i in 2:ns){
        points(density(phi[burnin:MCMC, iix$ix[i]]),
               type = "l",
               col = iix$ix[i],
               lwd = 2)
    }
    names <- c(sourcenames, rep("Unknown",
                                10))
    legend(x = "topright",
           names[iix$ix],
           text.col = iix$ix)
}
