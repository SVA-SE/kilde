plota <- function(mcmc_ob, burnin){
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
    
## dev.new()

## loASP<-matrix(0,mcmc_ob$var_a$ns,nat[1]); loGLN<-matrix(0,mcmc_ob$var_a$ns,nat[2])
## loGLT<-matrix(0,mcmc_ob$var_a$ns,nat[3]); loGLY<-matrix(0,mcmc_ob$var_a$ns,nat[4])
## loPGM<-matrix(0,mcmc_ob$var_a$ns,nat[5]); loTKT<-matrix(0,mcmc_ob$var_a$ns,nat[6])
## loUNC<-matrix(0,mcmc_ob$var_a$ns,nat[7])
## upASP<-matrix(0,mcmc_ob$var_a$ns,nat[1]); upGLN<-matrix(0,mcmc_ob$var_a$ns,nat[2])
## upGLT<-matrix(0,mcmc_ob$var_a$ns,nat[3]); upGLY<-matrix(0,mcmc_ob$var_a$ns,nat[4])
## upPGM<-matrix(0,mcmc_ob$var_a$ns,nat[5]); upTKT<-matrix(0,mcmc_ob$var_a$ns,nat[6])
## upUNC<-matrix(0,mcmc_ob$var_a$ns,nat[7])

## # check the fit of estimated freqs and observed freqs:
## # Set the number of picture frames as suitable:
## par(mfrow=c(2,4))
## errorASP <- matrix(NA,mcmc_ob$var_a$ns-1,nat[1]);errorGLN <- matrix(NA,mcmc_ob$var_a$ns-1,nat[2]);
## errorGLT <- matrix(NA,mcmc_ob$var_a$ns-1,nat[3]);errorGLY <- matrix(NA,mcmc_ob$var_a$ns-1,nat[4]);
## errorPGM <- matrix(NA,mcmc_ob$var_a$ns-1,nat[5]);errorTKT <- matrix(NA,mcmc_ob$var_a$ns-1,nat[6]);
## errorUNC <- matrix(NA,mcmc_ob$var_a$ns-1,nat[7]);
## errorASPhum <- numeric();errorGLNhum <- numeric();
## errorGLThum <- numeric();errorGLYhum <- numeric();
## errorPGMhum <- numeric();errorTKThum <- numeric();
## errorUNChum <- numeric();

## mqASP <- matrix(NA,mcmc_ob$var_a$ns,nat[1])
## for(i in 1:mcmc_ob$var_a$ns){for(j in 1:nat[1]){
## mqASP[i,j] <- mean(qASP[burnin:MCMC,i,j]) 
## loup <- quantile(qASP[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
## loASP[i,j] <- loup[1]; upASP[i,j] <- loup[2]
## } }
## plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="ASP")
## for(i in 1:(mcmc_ob$var_a$ns-1)){
## errorASP[i,1:nat[1]] <- ( mqASP[i,]-sourcesASP[i,]/sum(sourcesASP[i,]) )^2
## points(mqASP[i,],sourcesASP[i,]/sum(sourcesASP[i,]),col=i,pch=16)
## for(j in 1:nat[1]){points(c(loASP[i,j],upASP[i,j]),rep(sourcesASP[i,j]/sum(sourcesASP[i,]),2),col=i,type="l")}
## }
## mqGLN <- matrix(NA,mcmc_ob$var_a$ns,nat[2])
## for(i in 1:mcmc_ob$var_a$ns){for(j in 1:nat[2]){
## mqGLN[i,j] <- mean(qGLN[burnin:MCMC,i,j]) 
## loup <- quantile(qGLN[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
## loGLN[i,j] <- loup[1]; upGLN[i,j] <- loup[2]
## } }
## plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="GLN")
## for(i in 1:(mcmc_ob$var_a$ns-1)){
## errorGLN[i,1:nat[2]] <- ( mqGLN[i,]-sourcesGLN[i,]/sum(sourcesGLN[i,]) )^2
## points(mqGLN[i,],sourcesGLN[i,]/sum(sourcesGLN[i,]),col=i,pch=16)
## for(j in 1:nat[2]){points(c(loGLN[i,j],upGLN[i,j]),rep(sourcesGLN[i,j]/sum(sourcesGLN[i,]),2),col=i,type="l")}
## }
## mqGLT <- matrix(NA,mcmc_ob$var_a$ns,nat[3])
## for(i in 1:mcmc_ob$var_a$ns){for(j in 1:nat[3]){
## mqGLT[i,j] <- mean(qGLT[burnin:MCMC,i,j]) 
## loup <- quantile(qGLT[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
## loGLT[i,j] <- loup[1]; upGLT[i,j] <- loup[2]
## } }
## plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="GLT")
## for(i in 1:(mcmc_ob$var_a$ns-1)){
## errorGLT[i,1:nat[3]] <- ( mqGLT[i,]-sourcesGLT[i,]/sum(sourcesGLT[i,]) )^2
## points(mqGLT[i,],sourcesGLT[i,]/sum(sourcesGLT[i,]),col=i,pch=16)
## for(j in 1:nat[3]){points(c(loGLT[i,j],upGLT[i,j]),rep(sourcesGLT[i,j]/sum(sourcesGLT[i,]),2),col=i,type="l")}
## }
## mqGLY <- matrix(NA,mcmc_ob$var_a$ns,nat[4])
## for(i in 1:mcmc_ob$var_a$ns){for(j in 1:nat[4]){
## mqGLY[i,j] <- mean(qGLY[burnin:MCMC,i,j]) 
## loup <- quantile(qGLY[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
## loGLY[i,j] <- loup[1]; upGLY[i,j] <- loup[2]
## } }
## plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="GLY")
## for(i in 1:(mcmc_ob$var_a$ns-1)){
## errorGLY[i,1:nat[4]] <- ( mqGLY[i,]-sourcesGLY[i,]/sum(sourcesGLY[i,]) )^2
## points(mqGLY[i,],sourcesGLY[i,]/sum(sourcesGLY[i,]),col=i,pch=16)
## for(j in 1:nat[4]){points(c(loGLY[i,j],upGLY[i,j]),rep(sourcesGLY[i,j]/sum(sourcesGLY[i,]),2),col=i,type="l")}
## }
## mqPGM <- matrix(NA,mcmc_ob$var_a$ns,nat[5])
## for(i in 1:mcmc_ob$var_a$ns){for(j in 1:nat[5]){
## mqPGM[i,j] <- mean(qPGM[burnin:MCMC,i,j]) 
## loup <- quantile(qPGM[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
## loPGM[i,j] <- loup[1]; upPGM[i,j] <- loup[2]
## } }
## plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="PGM")
## for(i in 1:(mcmc_ob$var_a$ns-1)){
## errorPGM[i,1:nat[5]] <- ( mqPGM[i,]-sourcesPGM[i,]/sum(sourcesPGM[i,]) )^2
## points(mqPGM[i,],sourcesPGM[i,]/sum(sourcesPGM[i,]),col=i,pch=16)
## for(j in 1:nat[5]){points(c(loPGM[i,j],upPGM[i,j]),rep(sourcesPGM[i,j]/sum(sourcesPGM[i,]),2),col=i,type="l")}
## }
## mqTKT <- matrix(NA,mcmc_ob$var_a$ns,nat[6])
## for(i in 1:mcmc_ob$var_a$ns){for(j in 1:nat[6]){
## mqTKT[i,j] <- mean(qTKT[burnin:MCMC,i,j]) 
## loup <- quantile(qTKT[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
## loTKT[i,j] <- loup[1]; upTKT[i,j] <- loup[2]
## } }
## plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="TKT")
## for(i in 1:(mcmc_ob$var_a$ns-1)){
## errorTKT[i,1:nat[6]] <- ( mqTKT[i,]-sourcesTKT[i,]/sum(sourcesTKT[i,]) )^2
## points(mqTKT[i,],sourcesTKT[i,]/sum(sourcesTKT[i,]),col=i,pch=16)
## for(j in 1:nat[6]){points(c(loTKT[i,j],upTKT[i,j]),rep(sourcesTKT[i,j]/sum(sourcesTKT[i,]),2),col=i,type="l")}
## }
## mqUNC <- matrix(NA,mcmc_ob$var_a$ns,nat[7])
## for(i in 1:mcmc_ob$var_a$ns){for(j in 1:nat[7]){
## mqUNC[i,j] <- mean(qUNC[burnin:MCMC,i,j]) 
## loup <- quantile(qUNC[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
## loUNC[i,j] <- loup[1]; upUNC[i,j] <- loup[2]
## } }
## plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="UNC")
## for(i in 1:(mcmc_ob$var_a$ns-1)){
## errorUNC[i,1:nat[7]] <- ( mqUNC[i,]-sourcesUNC[i,]/sum(sourcesUNC[i,]) )^2
## points(mqUNC[i,],sourcesUNC[i,]/sum(sourcesUNC[i,]),col=i,pch=16)
## for(j in 1:nat[7]){points(c(loUNC[i,j],upUNC[i,j]),rep(sourcesUNC[i,j]/sum(sourcesUNC[i,]),2),col=i,type="l")}
## }

## # start plotting SA results:
## s <- numeric()
## for(i in 1:mcmc_ob$var_a$ns){ s[i] <- sd(mcmc_ob$var_a$phi[burnin:MCMC,i]) }
## iix <- sort(s,index.return=TRUE)
## attach(iix)
##   #pdf("SAresult.pdf")
## plot(density(mcmc_ob$var_a$phi[burnin:MCMC,ix[1]]),col=ix[1],xlim=c(0,1),xlab="Source attribution",ylab="",main=Ctr,lwd=2)
## for(i in 2:mcmc_ob$var_a$ns){
## points(density(mcmc_ob$var_a$phi[burnin:MCMC,ix[i]]),type="l",col=ix[i],lwd=2)
## }
## names <- c(sourcenames,rep("Unknown",10))
## legend(x="topright",names[ix],text.col=ix)
##   #dev.off()
## dev.new()
## hist(Z[burnin:MCMC,],0.5:(mcmc_ob$var_a$ns+0.5),col="turquoise3",labels=c(sourcenames,"Other"),axes=FALSE,ylab="",xlab="",main="P(Z | data)")


## errorsum <- numeric() 
## for(i in 1:mcmc_ob$var_a$ns-1){
## errorsum[i] <- sum(errorASP[i,])+sum(errorGLN[i,])+sum(errorGLT[i,])+sum(errorGLY[i,])+sum(errorPGM[i,])+sum(errorTKT[i,])+sum(errorUNC[i,])
## }
## mphi <- numeric()
## meanASPhum <- numeric();meanGLNhum <- numeric();meanGLThum <- numeric();meanGLYhum <- numeric()
## meanPGMhum <- numeric();meanTKThum <- numeric();meanUNChum <- numeric();

## for(i in 1:mcmc_ob$var_a$ns){mphi[i] <- mean(mcmc_ob$var_a$phi[burnin:MCMC,i]) } 
## for(j in 1:nat[1]){ meanASPhum[j] <- sum(mphi*mqASP[,j]);errorASPhum[j]<-(humansASP[j]/sum(humansASP)-meanASPhum[j])^2 }
## for(j in 1:nat[2]){ meanGLNhum[j] <- sum(mphi*mqGLN[,j]);errorGLNhum[j]<-(humansGLN[j]/sum(humansGLN)-meanGLNhum[j])^2 }
## for(j in 1:nat[3]){ meanGLThum[j] <- sum(mphi*mqGLT[,j]);errorGLThum[j]<-(humansGLT[j]/sum(humansGLT)-meanGLThum[j])^2 }
## for(j in 1:nat[4]){ meanGLYhum[j] <- sum(mphi*mqGLY[,j]);errorGLYhum[j]<-(humansGLY[j]/sum(humansGLY)-meanGLYhum[j])^2 }
## for(j in 1:nat[5]){ meanPGMhum[j] <- sum(mphi*mqPGM[,j]);errorPGMhum[j]<-(humansPGM[j]/sum(humansPGM)-meanPGMhum[j])^2 }
## for(j in 1:nat[6]){ meanTKThum[j] <- sum(mphi*mqTKT[,j]);errorTKThum[j]<-(humansTKT[j]/sum(humansTKT)-meanTKThum[j])^2 }
## for(j in 1:nat[7]){ meanUNChum[j] <- sum(mphi*mqUNC[,j]);errorUNChum[j]<-(humansUNC[j]/sum(humansUNC)-meanUNChum[j])^2 }

## ES <- sum(errorsum)   # sum of all squared errors
## EShum <- sum(errorASPhum)+sum(errorGLNhum)+sum(errorGLThum)+sum(errorGLYhum)+sum(errorPGMhum)+sum(errorTKThum)+sum(errorUNChum)


