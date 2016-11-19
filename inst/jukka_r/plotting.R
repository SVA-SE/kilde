plot_r_mcmc <- function(mcmc_ob, burnin){
    dev.new()
    par(mfrow = c(mcmc_ob$var_a$ns, 1), mar = c(2.5, 5.5, 2.5, 2.5))
    for(i in 1 : mcmc_ob$var_a$ns){
        plot(mcmc_ob$var_a$phi[burnin:mcmc_ob$var_a$MCMC, i],
             xlab="",
             ylab=mcmc_ob$var_b$data$sourcenames[i],
             pch=16,
             cex=0.5,
             cex.lab=1.7
             )
    }
}
    
loASP <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[1])
loGLN <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[2])
loGLT <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[3])
loGLY <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[4])
loPGM <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[5])
loTKT <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[6])
loUNC <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[7])
##
upASP <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[1])
upGLN <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[2])
upGLT <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[3])
upGLY <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[4])
upPGM <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[5])
upTKT <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[6])
upUNC <- matrix(0, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[7])
##
# check the fit of estimated freqs and observed freqs:
# Set the number of picture frames as suitable:
par(mfrow=c(2,4))
errorASP <- matrix(NA, mcmc_ob$var_a$ns - 1, mcmc_ob$var_b$inits$nat[1])
errorGLN <- matrix(NA, mcmc_ob$var_a$ns - 1, mcmc_ob$var_b$inits$nat[2])
errorGLT <- matrix(NA, mcmc_ob$var_a$ns - 1, mcmc_ob$var_b$inits$nat[3])
errorGLY <- matrix(NA, mcmc_ob$var_a$ns - 1, mcmc_ob$var_b$inits$nat[4])
errorPGM <- matrix(NA, mcmc_ob$var_a$ns - 1, mcmc_ob$var_b$inits$nat[5])
errorTKT <- matrix(NA, mcmc_ob$var_a$ns - 1, mcmc_ob$var_b$inits$nat[6])
errorUNC <- matrix(NA, mcmc_ob$var_a$ns - 1, mcmc_ob$var_b$inits$nat[7])
errorASPhum <- numeric()
errorGLNhum <- numeric()
errorGLThum <- numeric()
errorGLYhum <- numeric()
errorPGMhum <- numeric()
errorTKThum <- numeric()
errorUNChum <- numeric()
## ASP plot
mqASP <- matrix(NA, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[1])
for(i in 1:mcmc_ob$var_a$ns){
    for(j in 1:mcmc_ob$var_b$inits$nat[1]){
        mqASP[i, j] <- mean(mcmc_ob$var_a$qASP[burnin:mcmc_ob$var_a$MCMC, i, j]) 
        loup <- quantile(mcmc_ob$var_a$qASP[burnin:mcmc_ob$var_a$MCMC, i, j], c(0.025, 0.975), names=FALSE)
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
for(i in 1:(mcmc_ob$var_a$ns - 1)){
    errorASP[i, 1:mcmc_ob$var_b$inits$nat[1]] <- (mqASP[i, ] - mcmc_ob$var_b$data$sourcesASP[i, ] /
                              sum(mcmc_ob$var_b$data$sourcesASP[i, ]))^2
    points(mqASP[i, ], mcmc_ob$var_b$data$sourcesASP[i, ] /
                       sum(mcmc_ob$var_b$data$sourcesASP[i, ]), col=i, pch=16)
    for(j in 1:mcmc_ob$var_b$inits$nat[1]){ 
        points(c(loASP[i, j], upASP[i, j]),
               rep(mcmc_ob$var_b$data$sourcesASP[i, j] /
                   sum(mcmc_ob$var_b$data$sourcesASP[i, ]), 2), col=i, type="l")
    }
}
## GLN plot
mqGLN <- matrix(NA, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[2])
for(i in 1:mcmc_ob$var_a$ns){
    for(j in 1:mcmc_ob$var_b$inits$nat[2]){
        mqGLN[i, j] <- mean(mcmc_ob$var_a$qGLN[burnin:mcmc_ob$var_a$MCMC, i, j]) 
        loup <- quantile(mcmc_ob$var_a$qGLN[burnin:mcmc_ob$var_a$MCMC, i, j], c(0.025, 0.975), names=FALSE)
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
for(i in 1:(mcmc_ob$var_a$ns - 1)){
    errorGLN[i, 1:mcmc_ob$var_b$inits$nat[2]] <- (mqGLN[i, ] - mcmc_ob$var_b$data$sourcesGLN[i, ] /
                              sum(mcmc_ob$var_b$data$sourcesGLN[i, ]))^2
    points(mqGLN[i, ], mcmc_ob$var_b$data$sourcesGLN[i, ] /
                       sum(mcmc_ob$var_b$data$sourcesGLN[i, ]), col=i, pch=16)
    for(j in 1:mcmc_ob$var_b$inits$nat[2]){ 
        points(c(loGLN[i, j], upGLN[i, j]),
               rep(mcmc_ob$var_b$data$sourcesGLN[i, j] /
                   sum(mcmc_ob$var_b$data$sourcesGLN[i, ]), 2), col=i, type="l")
    }
}
## GLT plot
mqGLT <- matrix(NA, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[3])
for(i in 1:mcmc_ob$var_a$ns){
    for(j in 1:mcmc_ob$var_b$inits$nat[3]){
        mqGLT[i, j] <- mean(mcmc_ob$var_a$qGLT[burnin:mcmc_ob$var_a$MCMC, i, j]) 
        loup <- quantile(mcmc_ob$var_a$qGLT[burnin:mcmc_ob$var_a$MCMC, i, j], c(0.025, 0.975), names=FALSE)
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
for(i in 1:(mcmc_ob$var_a$ns - 1)){
    errorGLT[i, 1:mcmc_ob$var_b$inits$nat[3]] <- (mqGLT[i, ] - mcmc_ob$var_b$data$sourcesGLT[i, ] /
                              sum(mcmc_ob$var_b$data$sourcesGLT[i, ]))^2
    points(mqGLT[i, ], mcmc_ob$var_b$data$sourcesGLT[i, ] /
                       sum(mcmc_ob$var_b$data$sourcesGLT[i, ]), col=i, pch=16)
    for(j in 1:mcmc_ob$var_b$inits$nat[3]){ 
        points(c(loGLT[i, j], upGLT[i, j]),
               rep(mcmc_ob$var_b$data$sourcesGLT[i, j] /
                   sum(mcmc_ob$var_b$data$sourcesGLT[i, ]), 2), col=i, type="l")
    }
}
## GLY plot
mqGLY <- matrix(NA, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[4])
for(i in 1:mcmc_ob$var_a$ns){
    for(j in 1:mcmc_ob$var_b$inits$nat[4]){
        mqGLY[i, j] <- mean(mcmc_ob$var_a$qGLY[burnin:mcmc_ob$var_a$MCMC, i, j]) 
        loup <- quantile(mcmc_ob$var_a$qGLY[burnin:mcmc_ob$var_a$MCMC, i, j], c(0.025, 0.975), names=FALSE)
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
for(i in 1:(mcmc_ob$var_a$ns - 1)){
    errorGLY[i, 1:mcmc_ob$var_b$inits$nat[4]] <- (mqGLY[i, ] - mcmc_ob$var_b$data$sourcesGLY[i, ] /
                              sum(mcmc_ob$var_b$data$sourcesGLY[i, ]))^2
    points(mqGLY[i, ], mcmc_ob$var_b$data$sourcesGLY[i, ] /
                       sum(mcmc_ob$var_b$data$sourcesGLY[i, ]), col=i, pch=16)
    for(j in 1:mcmc_ob$var_b$inits$nat[4]){ 
        points(c(loGLY[i, j], upGLY[i, j]),
               rep(mcmc_ob$var_b$data$sourcesGLY[i, j] /
                   sum(mcmc_ob$var_b$data$sourcesGLY[i, ]), 2), col=i, type="l")
    }
}
## PGM plot
mqPGM <- matrix(NA, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[5])
for(i in 1:mcmc_ob$var_a$ns){
    for(j in 1:mcmc_ob$var_b$inits$nat[5]){
        mqPGM[i, j] <- mean(mcmc_ob$var_a$qPGM[burnin:mcmc_ob$var_a$MCMC, i, j]) 
        loup <- quantile(mcmc_ob$var_a$qPGM[burnin:mcmc_ob$var_a$MCMC, i, j], c(0.025, 0.975), names=FALSE)
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
for(i in 1:(mcmc_ob$var_a$ns - 1)){
    errorPGM[i, 1:mcmc_ob$var_b$inits$nat[5]] <- (mqPGM[i, ] - mcmc_ob$var_b$data$sourcesPGM[i, ] /
                              sum(mcmc_ob$var_b$data$sourcesPGM[i, ]))^2
    points(mqPGM[i, ], mcmc_ob$var_b$data$sourcesPGM[i, ] /
                       sum(mcmc_ob$var_b$data$sourcesPGM[i, ]), col=i, pch=16)
    for(j in 1:mcmc_ob$var_b$inits$nat[5]){ 
        points(c(loPGM[i, j], upPGM[i, j]),
               rep(mcmc_ob$var_b$data$sourcesPGM[i, j] /
                   sum(mcmc_ob$var_b$data$sourcesPGM[i, ]), 2), col=i, type="l")
    }
}
## TKT plot
mqTKT <- matrix(NA, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[6])
for(i in 1:mcmc_ob$var_a$ns){
    for(j in 1:mcmc_ob$var_b$inits$nat[6]){
        mqTKT[i, j] <- mean(mcmc_ob$var_a$qTKT[burnin:mcmc_ob$var_a$MCMC, i, j]) 
        loup <- quantile(mcmc_ob$var_a$qTKT[burnin:mcmc_ob$var_a$MCMC, i, j], c(0.025, 0.975), names=FALSE)
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
for(i in 1:(mcmc_ob$var_a$ns - 1)){
    errorTKT[i, 1:mcmc_ob$var_b$inits$nat[6]] <- (mqTKT[i, ] - mcmc_ob$var_b$data$sourcesTKT[i, ] /
                              sum(mcmc_ob$var_b$data$sourcesTKT[i, ]))^2
    points(mqTKT[i, ], mcmc_ob$var_b$data$sourcesTKT[i, ] /
                       sum(mcmc_ob$var_b$data$sourcesTKT[i, ]), col=i, pch=16)
    for(j in 1:mcmc_ob$var_b$inits$nat[6]){ 
        points(c(loTKT[i, j], upTKT[i, j]),
               rep(mcmc_ob$var_b$data$sourcesTKT[i, j] /
                   sum(mcmc_ob$var_b$data$sourcesTKT[i, ]), 2), col=i, type="l")
    }
}
## UNC plot
mqUNC <- matrix(NA, mcmc_ob$var_a$ns, mcmc_ob$var_b$inits$nat[7])
for(i in 1:mcmc_ob$var_a$ns){
    for(j in 1:mcmc_ob$var_b$inits$nat[7]){
        mqUNC[i, j] <- mean(mcmc_ob$var_a$qUNC[burnin:mcmc_ob$var_a$MCMC, i, j]) 
        loup <- quantile(mcmc_ob$var_a$qUNC[burnin:mcmc_ob$var_a$MCMC, i, j], c(0.025, 0.975), names=FALSE)
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
for(i in 1:(mcmc_ob$var_a$ns - 1)){
    errorUNC[i, 1:mcmc_ob$var_b$inits$nat[7]] <- (mqUNC[i, ] - mcmc_ob$var_b$data$sourcesUNC[i, ] /
                              sum(mcmc_ob$var_b$data$sourcesUNC[i, ]))^2
    points(mqUNC[i, ], mcmc_ob$var_b$data$sourcesUNC[i, ] /
                       sum(mcmc_ob$var_b$data$sourcesUNC[i, ]), col=i, pch=16)
    for(j in 1:mcmc_ob$var_b$inits$nat[7]){ 
        points(c(loUNC[i, j], upUNC[i, j]),
               rep(mcmc_ob$var_b$data$sourcesUNC[i, j] /
                   sum(mcmc_ob$var_b$data$sourcesUNC[i, ]), 2), col=i, type="l")
    }
}
#####



# start plotting SA results:
s <- numeric()
for(i in 1:mcmc_ob$var_a$ns){
    s[i] <- sd(mcmc_ob$var_a$phi[burnin:mcmc_ob$var_a$MCMC, i])
}
iix <- sort(s, index.return = TRUE)
attach(iix)
##pdf("SAresult.pdf")
plot(density(mcmc_ob$var_a$phi[burnin:mcmc_ob$var_a$MCMC, iix$ix[1]]),
     col = iix$ix[1],
     xlim = c(0, 1),
     xlab = "Source attribution",
     ylab = "",
     main = "",
     lwd = 2)
for(i in 2:mcmc_ob$var_a$ns){
    points(density(mcmc_ob$var_a$phi[burnin:mcmc_ob$var_a$MCMC, iix$ix[i]]),
           type = "l",
           col = iix$ix[i],
           lwd = 2)
}
names <- c(mcmc_ob$var_b$data$sourcenames, rep("Unknown", 10))
legend(x = "topright",
       names[iix$ix],
       text.col=iix$ix)

dev.new()
hist(Z[burnin:MCMC, ],
     0.5:(mcmc_ob$var_a$ns + 0.5),
     col = "turquoise3",
     labels = c(sourcenames, "Other"),
     axes = FALSE,
     ylab = "",
     xlab = "",
     main = "P(Z | data)")
#############
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


