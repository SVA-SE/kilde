
function(phi, ns, MCMC, burnin, sourcenames){
    savepar <- par(mfrow = c(ns, 1),
                   mar = c(2.5, 5.5, 2.5, 2.5)
                   )
    on.exit(par(savepar))
    for(i in 1:ns){
        plot(phi[burnin:MCMC, i],
             xlab = "",
             ylab = sourcenames[i],
             pch = 16,
             cex = 0.4,
             cex.lab = 1.7)
    }
}


dev.new()

loST<-matrix(0,ns,length(STu))
upST<-matrix(0,ns,length(STu))

# check the fit of estimated freqs and observed freqs:
# Set the number of picture frames as suitable:
par(mfrow=c(2,1))
errorST <- matrix(NA,ns-1,length(STu))
errorSThum <- numeric()
mqST <- matrix(NA,ns,length(STu))
fST <- numeric()
for(i in 1:(ns-1)){for(j in 1:length(STu)){fST[j] <- sourcesST[i,j]/sum(sourcesST[i,]) }}
for(i in 1:ns){for(j in 1:length(STu)){
mqST[i,j] <- mean(qST[burnin:MCMC,i,j]) 
loup <- quantile(qST[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
loST[i,j] <- loup[1]; upST[i,j] <- loup[2]
} }
limx <- max(upST)*1.01
limy <- max(fST)*1.1 
plot(c(0,1),c(0,1),type="l",xlim=c(0,limx),ylim=c(0,limy),xlab="Estimated freq",ylab="Observed sample freq",main="ST")
for(i in 1:(ns-1)){
errorST[i,1:length(STu)] <- ( mqST[i,]-sourcesST[i,]/sum(sourcesST[i,]) )^2
points(mqST[i,],sourcesST[i,]/sum(sourcesST[i,]),col=i,pch=16)
for(j in 1:length(STu)){points(c(loST[i,j],upST[i,j]),rep(sourcesST[i,j]/sum(sourcesST[i,]),2),col=i,type="l")}
}
# start plotting SA results:
s <- numeric()
for(i in 1:ns){ s[i] <- sd(phi[burnin:MCMC,i]) }
iix <- sort(s,index.return=TRUE)
attach(iix)
  #pdf("SAresult.pdf")
plot(density(phi[burnin:MCMC,ix[1]]),ylim=c(0,60),col=ix[1],xlim=c(0,1),xlab="Source attribution",ylab="",main=Ctr,lwd=2)
for(i in 2:ns){
points(density(phi[burnin:MCMC,ix[i]]),type="l",col=ix[i],lwd=2)
}
names <- c(sourcenames,rep("Unknown",10))
legend(x="topright",names[ix],text.col=ix)
  #dev.off()
dev.new()
hist(Z[burnin:MCMC,],0.5:(ns+0.5),col="turquoise3",labels=c(sourcenames,"Other"),axes=FALSE,ylab="",xlab="",main="P(Z | data)")

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

