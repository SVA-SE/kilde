##' summary_internal
##'
##' A summary of a model run
##' @title summary_internal
##' @param ns ns 
##' @param humansASP humansASP 
##' @param humansGLN humansGLN 
##' @param humansGLT humansGLT 
##' @param humansGLY humansGLY 
##' @param humansPGM humansPGM 
##' @param humansTKT humansTKT 
##' @param humansUNC humansUNC 
##' @param sourcesASP sourcesASP 
##' @param sourcesGLN sourcesGLN 
##' @param sourcesGLT sourcesGLT 
##' @param sourcesGLY sourcesGLY 
##' @param sourcesPGM sourcesPGM 
##' @param sourcesTKT sourcesTKT 
##' @param sourcesUNC sourcesUNC 
##' @param qASP qASP
##' @param qGLN qGLN
##' @param qGLT qGLT
##' @param qGLY qGLY
##' @param qPGM qPGM
##' @param qTKT qTKT
##' @param qUNC qUNC
##' @param phi phi 
##' @param burnin burnin 
##' @param MCMC MCMC 
##' @param nat nat 
##' @return A summary
##' @author Jukka Ranta
summary_internal <- function(ns,
                             humansASP,
                             humansGLN,
                             humansGLT,
                             humansGLY,
                             humansPGM,
                             humansTKT,
                             humansUNC,
                             sourcesASP,
                             sourcesGLN,
                             sourcesGLT,
                             sourcesGLY,
                             sourcesPGM,
                             sourcesTKT,
                             sourcesUNC,
                             qASP,
                             qGLN,
                             qGLT,
                             qGLY,
                             qPGM,
                             qTKT,
                             qUNC,
                             phi,
                             burnin,
                             MCMC,
                             nat
                             ){
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
    ## ASP summary
    mqASP <- matrix(NA, ns, nat[1])
    for(i in 1:ns){
        for(j in 1:nat[1]){
            mqASP[i, j] <- mean(qASP[burnin:MCMC, i, j]) 
            loup <- quantile(qASP[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loASP[i, j] <- loup[1]; upASP[i,j] <- loup[2]
        }
    }
    for(i in 1:(ns - 1)){
        errorASP[i, 1:nat[1]] <- (mqASP[i, ] - sourcesASP[i, ] /
                                  sum(sourcesASP[i, ]))^2
    }
    ## GLN summary
    mqGLN <- matrix(NA, ns, nat[2])
    for(i in 1:ns){
        for(j in 1:nat[2]){
            mqGLN[i, j] <- mean(qGLN[burnin:MCMC, i, j]) 
            loup <- quantile(qGLN[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loGLN[i, j] <- loup[1]; upGLN[i,j] <- loup[2]
        }
    }
    for(i in 1:(ns - 1)){
        errorGLN[i, 1:nat[2]] <- (mqGLN[i, ] - sourcesGLN[i, ] /
                                  sum(sourcesGLN[i, ]))^2
    }
    ## GLT summary
    mqGLT <- matrix(NA, ns, nat[3])
    for(i in 1:ns){
        for(j in 1:nat[3]){
            mqGLT[i, j] <- mean(qGLT[burnin:MCMC, i, j]) 
            loup <- quantile(qGLT[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loGLT[i, j] <- loup[1]; upGLT[i,j] <- loup[2]
        }
    }
    for(i in 1:(ns - 1)){
        errorGLT[i, 1:nat[3]] <- (mqGLT[i, ] - sourcesGLT[i, ] /
                                  sum(sourcesGLT[i, ]))^2
    }
    ## GLY summary
    mqGLY <- matrix(NA, ns, nat[4])
    for(i in 1:ns){
        for(j in 1:nat[4]){
            mqGLY[i, j] <- mean(qGLY[burnin:MCMC, i, j]) 
            loup <- quantile(qGLY[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loGLY[i, j] <- loup[1]; upGLY[i,j] <- loup[2]
        }
    }
    for(i in 1:(ns - 1)){
        errorGLY[i, 1:nat[4]] <- (mqGLY[i, ] - sourcesGLY[i, ] /
                                  sum(sourcesGLY[i, ]))^2
    }
    ## PGM summary
    mqPGM <- matrix(NA, ns, nat[5])
    for(i in 1:ns){
        for(j in 1:nat[5]){
            mqPGM[i, j] <- mean(qPGM[burnin:MCMC, i, j]) 
            loup <- quantile(qPGM[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loPGM[i, j] <- loup[1]; upPGM[i,j] <- loup[2]
        }
    }
    for(i in 1:(ns - 1)){
        errorPGM[i, 1:nat[5]] <- (mqPGM[i, ] - sourcesPGM[i, ] /
                                                      sum(sourcesPGM[i, ]))^2
    }
    ## TKT summary
    mqTKT <- matrix(NA, ns, nat[6])
    for(i in 1:ns){
        for(j in 1:nat[6]){
            mqTKT[i, j] <- mean(qTKT[burnin:MCMC, i, j]) 
            loup <- quantile(qTKT[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loTKT[i, j] <- loup[1]; upTKT[i,j] <- loup[2]
        }
    }
    for(i in 1:(ns - 1)){
        errorTKT[i, 1:nat[6]] <- (mqTKT[i, ] - sourcesTKT[i, ] /
                                                      sum(sourcesTKT[i, ]))^2
    }
    ## UNC summary
    mqUNC <- matrix(NA, ns, nat[7])
    for(i in 1:ns){
        for(j in 1:nat[7]){
            mqUNC[i, j] <- mean(qUNC[burnin:MCMC, i, j]) 
            loup <- quantile(qUNC[burnin:MCMC, i, j], c(0.025, 0.975), names=FALSE)
            loUNC[i, j] <- loup[1]; upUNC[i,j] <- loup[2]
        }
    }
    for(i in 1:(ns - 1)){
        errorUNC[i, 1:nat[7]] <- (mqUNC[i, ] - sourcesUNC[i, ] /
                                                      sum(sourcesUNC[i, ]))^2
    }
    ## start plotting SA results:
    s <- numeric()
    for(i in 1:ns){
        s[i] <- sd(phi[burnin:MCMC, i])
    }
    iix <- sort(s, index.return = TRUE)
################################
################################
    errorsum <- numeric() 
    for(i in 1:ns-1){
        errorsum[i] <-
            sum(errorASP[i, ]) +
            sum(errorGLN[i, ]) +
            sum(errorGLT[i, ]) +
            sum(errorGLY[i, ]) +
            sum(errorPGM[i, ]) +
            sum(errorTKT[i, ]) +
            sum(errorUNC[i, ])
    }
    mphi <- numeric()
    meanASPhum <- numeric()
    meanGLNhum <- numeric()
    meanGLThum <- numeric()
    meanGLYhum <- numeric()
    meanPGMhum <- numeric()
    meanTKThum <- numeric()
    meanUNChum <- numeric()
    for(i in 1:ns){
        mphi[i] <- mean(phi[burnin:MCMC, i])
    } 
    for(j in 1:nat[1]){
        meanASPhum[j] <- sum(mphi * mqASP[, j])
        errorASPhum[j] <- (humansASP[j] / sum(humansASP) - meanASPhum[j])^2
    }
    for(j in 1:nat[2]){
        meanGLNhum[j] <- sum(mphi * mqGLN[, j])
        errorGLNhum[j] <- (humansGLN[j] / sum(humansGLN) - meanGLNhum[j])^2
    }
    for(j in 1:nat[3]){
        meanGLThum[j] <- sum(mphi * mqGLT[, j])
        errorGLThum[j] <- (humansGLT[j] / sum(humansGLT) - meanGLThum[j])^2
    }
    for(j in 1:nat[4]){
        meanGLYhum[j] <- sum(mphi * mqGLY[, j])
        errorGLYhum[j] <- (humansGLY[j] / sum(humansGLY) - meanGLYhum[j])^2
    }
    for(j in 1:nat[5]){
        meanPGMhum[j] <- sum(mphi * mqPGM[, j])
        errorPGMhum[j] <- (humansPGM[j] / sum(humansPGM) - meanPGMhum[j])^2
    }
    for(j in 1:nat[6]){
        meanTKThum[j] <- sum(mphi * mqTKT[, j])
        errorTKThum[j] <- (humansTKT[j] / sum(humansTKT) - meanTKThum[j])^2
    }
    for(j in 1:nat[7]){
        meanUNChum[j] <- sum(mphi * mqUNC[, j])
        errorUNChum[j] <- (humansUNC[j] / sum(humansUNC) - meanUNChum[j])^2
    }

    ES <- sum(errorsum)   # sum of all squared errors
    EShum <-
        sum(errorASPhum) +
        sum(errorGLNhum) +
        sum(errorGLThum) +
        sum(errorGLYhum) +
        sum(errorPGMhum) +
        sum(errorTKThum) +
        sum(errorUNChum)
    result <- list(ES = ES,
                   EShum = EShum)
    return(result)
}
##' summary of models
##' 
##' Defining the generic function
##' @title summary 
##' @param x The object to be summarized
##' @param burnin the burnin length
##' @return A list
##' @export
##' @author Thomas Rosendal
summary_kilde <- function(x, burnin){
    UseMethod('summary_kilde')
}
summary_kilde.default = summary_kilde
##' A summary method for the kilde_rmcmc object
##'
##' @title summary_kilde.kilde_rmcmc
##' @param x kilde_rmcmc object
##' @param burnin The burnin length
##' @return A summary
##' @export
##' @author Thomas Rosendal
summary_kilde.kilde_rmcmc <- function(x, burnin){
    ns <- x$var_b$inits$ns
    humansASP <- x$var_b$data$humansASP
    humansGLN <- x$var_b$data$humansGLN
    humansGLT <- x$var_b$data$humansGLT
    humansGLY <- x$var_b$data$humansGLY
    humansPGM <- x$var_b$data$humansPGM
    humansTKT <- x$var_b$data$humansTKT
    humansUNC <- x$var_b$data$humansUNC
    sourcesASP <- x$var_b$data$sourcesASP
    sourcesGLN <- x$var_b$data$sourcesGLN
    sourcesGLT <- x$var_b$data$sourcesGLT
    sourcesGLY <- x$var_b$data$sourcesGLY
    sourcesPGM <- x$var_b$data$sourcesPGM
    sourcesTKT <- x$var_b$data$sourcesTKT
    sourcesUNC <- x$var_b$data$sourcesUNC
    qASP <- x$var_a$qASP  
    qGLN <- x$var_a$qGLN  
    qGLT <- x$var_a$qGLT  
    qGLY <- x$var_a$qGLY  
    qPGM <- x$var_a$qPGM  
    qTKT <- x$var_a$qTKT  
    qUNC <- x$var_a$qUNC
    phi <- x$var_a$phi
    MCMC <- x$var_a$MCMC
    nat <- x$var_a$nat
    summary_internal(ns,
                     humansASP,
                     humansGLN,
                     humansGLT,
                     humansGLY,
                     humansPGM,
                     humansTKT,
                     humansUNC,
                     sourcesASP,
                     sourcesGLN,
                     sourcesGLT,
                     sourcesGLY,
                     sourcesPGM,
                     sourcesTKT,
                     sourcesUNC,
                     qASP,
                     qGLN,
                     qGLT,
                     qGLY,
                     qPGM,
                     qTKT,
                     qUNC,
                     phi,
                     burnin,
                     MCMC,
                     nat
                     )
}
##' A summary method for the bugsmcmc function
##'
##' @title summary_kilde.kilde_bugsmcmc
##' @param x A kilde_bugsmcmc object
##' @param burnin the burnin length
##' @return A summary
##' @export
##' @author Thomas Rosendal
summary_kilde.kilde_bugsmcmc <- function(x, burnin){
    ns <- x$var_b$inits$ns
    humansASP <- x$var_b$data$humansASP
    humansGLN <- x$var_b$data$humansGLN
    humansGLT <- x$var_b$data$humansGLT
    humansGLY <- x$var_b$data$humansGLY
    humansPGM <- x$var_b$data$humansPGM
    humansTKT <- x$var_b$data$humansTKT
    humansUNC <- x$var_b$data$humansUNC
    sourcesASP <- x$var_b$data$sourcesASP
    sourcesGLN <- x$var_b$data$sourcesGLN
    sourcesGLT <- x$var_b$data$sourcesGLT
    sourcesGLY <- x$var_b$data$sourcesGLY
    sourcesPGM <- x$var_b$data$sourcesPGM
    sourcesTKT <- x$var_b$data$sourcesTKT
    sourcesUNC <- x$var_b$data$sourcesUNC
    qASP <- x$bugs_result$sims.list$qASP  
    qGLN <- x$bugs_result$sims.list$qGLN  
    qGLT <- x$bugs_result$sims.list$qGLT  
    qGLY <- x$bugs_result$sims.list$qGLY  
    qPGM <- x$bugs_result$sims.list$qPGM  
    qTKT <- x$bugs_result$sims.list$qTKT  
    qUNC <- x$bugs_result$sims.list$qUNC
    phi <- x$bugs_result$sims.list$phi
    MCMC <- x$bugs_result$n.iter
    nat <- x$var_b$inits$nat
    MCMC <- MCMC - burnin
    burnin <- 1
    summary_internal(ns,
                     humansASP,
                     humansGLN,
                     humansGLT,
                     humansGLY,
                     humansPGM,
                     humansTKT,
                     humansUNC,
                     sourcesASP,
                     sourcesGLN,
                     sourcesGLT,
                     sourcesGLY,
                     sourcesPGM,
                     sourcesTKT,
                     sourcesUNC,
                     qASP,
                     qGLN,
                     qGLT,
                     qGLY,
                     qPGM,
                     qTKT,
                     qUNC,
                     phi,
                     burnin,
                     MCMC,
                     nat
                     )
}
##' summary_kilde.kilde_bugsmcmc_ST
##'
##' @title summary_kilde.kilde_bugsmcmc_ST 
##' @param x An object of class kilde_bugsmcmc_ST
##' @param burnin the burnin length 
##' @return A list
##' @export
##' @author Thomas Rosendal
summary_kilde.kilde_bugsmcmc_ST <- function(x, burnin){
    ns <- x$other$ns
    phi <- x$bugs_result$sims.list$phi
    MCMC <- x$bugs_result$n.iter
    humansST <- x$var_b$humansST
    STu <- x$var_b$STu
    sourcesST <- x$var_b$sourcesST
    qST <- x$bugs_result$sims.list$qST
    MCMC <- MCMC - burnin
    burnin <- 0
    summary_kilde_ST_internal(ns,
                              phi,
                              burnin,
                              MCMC,
                              humansST,
                              STu,
                              sourcesST,
                              qST)
}
##' summary_kilde.kilde_rmcmc_ST
##'
##' @title summary_kilde.kilde_rmcmc_ST 
##' @param x An object of class kilde_rmcmc_ST
##' @param burnin the burnin length 
##' @return A list
##' @export
##' @author Thomas Rosendal
summary_kilde.kilde_rmcmc_ST <- function(x, burnin){
    ns <- x$var_a$ns
    phi <- x$var_a$phi
    MCMC <- x$var_a$MCMC
    humansST <- x$var_b$humansST
    STu <- x$var_a$STu
    sourcesST <- x$var_b$sourcesST
    qST <- x$var_a$qST
    summary_kilde_ST_internal(ns,
                              phi,
                              burnin,
                              MCMC,
                              humansST,
                              STu,
                              sourcesST,
                              qST)
}
##' summary_kilde_ST_internal
##'
##' @title summary_kilde_ST_internal 
##' @param ns ns 
##' @param phi phi 
##' @param burnin burnin 
##' @param MCMC MCMC 
##' @param humansST humansST 
##' @param STu STu 
##' @param sourcesST sourcesST 
##' @param qST qST 
##' @return A list
##' @author Thomas Rosendal
summary_kilde_ST_internal <- function(ns,
                                      phi,
                                      burnin,
                                      MCMC,
                                      humansST,
                                      STu,
                                      sourcesST,
                                      qST){
    ## The population attribution plot and the fit
    loST<-matrix(0, ns, length(STu))
    upST<-matrix(0, ns, length(STu))
    ## check the fit of estimated freqs and observed freqs:
    ## Set the number of picture frames as suitable:
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
    for(i in 1:(ns - 1)){
        errorST[i, 1:length(STu)] <- ( mqST[i, ] - sourcesST[i, ] / sum(sourcesST[i, ]))^2
    }
    s <- numeric()
    for(i in 1:ns){
        s[i] <- sd(phi[burnin:MCMC, i])
    }
    iix <- sort(s, index.return=TRUE)
    errorsum <- numeric() 
    for(i in 1:ns - 1){
        errorsum[i] <- sum(errorST[i, ])
    }
    mphi <- numeric()
    meanSThum <- numeric()
    for(i in 1:ns){
        mphi[i] <- mean(phi[burnin:MCMC, i])
    } 
    for(j in 1:length(STu)){
        meanSThum[j] <- sum(mphi * mqST[, j])
        errorSThum[j] <- (humansST[j] / sum(humansST) - meanSThum[j])^2
    }   
    ES <- sum(errorsum)   # sum of all squared errors
    EShum <- sum(errorSThum)
    
    result <- list(ES = ES,
                   EShum = EShum)
    return(result)
}
