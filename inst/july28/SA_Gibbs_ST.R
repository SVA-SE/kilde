
DATA <- read.table("NMDD2015data_vers21_modif.txt",header=TRUE)

attach(DATA)
# use ALL allele types
z <- !is.na(ST)

Ctr = "SE"    # choose the Country
sourcenames=setdiff(unique(Group[Country==Ctr]),"Human")
#sourcenames=c("Turkey","Cattle","Broiler","Water") # FI  
#sourcenames=c("Bird","Cattle","Sheep","Broiler","Pig","Dog","Water") # SE 
#sourcenames=c("Cattle","Broiler","Pig","Turkey","Duck") # DK
#sourcenames=c("Broiler","Sheep","Pig","Cattle") # NO
ns <- length(sourcenames)

# find how many different STs there are (in all isolates), and list them:
STu <- sort(unique(ST[(Country==Ctr)&(z)]))
STuH <- sort(unique(ST[(Country==Ctr)&(z)&(Group=="Human")]))
STuS <- sort(unique(ST[(Country==Ctr)&(z)&(Group!="Human")]))
STuHo <- setdiff(STuH,STuS) 

HumanST <- ST[(z)&(Group=="Human")&(Country==Ctr)]
Humnovel <- 0
for(i in 1:length(HumanST)){
Humnovel <- Humnovel + is.element(HumanST[i],STuHo)*1
}
# Sample probability to see ST in Human isolates that was not seen in sources: 
PUNST <- Humnovel/length(HumanST)


# number (counts) of ST-types in each source:
# define "unknown source" (ns+1) with zero observed counts.
# OR: define many unknown sources: ns+K for the DPP model?
maxns <- ns+1
sourcesST  <- matrix(0,maxns,length(STu)) 

# number (counts) of such ST-types in humans & also found in sources:
humansST <- numeric()

# sample frequencies (counts) of each ST-type
for(i in 1:length(STu)){
humansST[i] <- sum((ST==STu[i])&(Group=="Human")&(Country==Ctr)&z)
for(j in 1:ns){
sourcesST[j,i] <- sum((ST==STu[i])&(Group==sourcenames[j])&(Country==Ctr)&z) }
}

# Number of all human isolates:
Nisolates <- sum((Country==Ctr)&(z)&(Group=="Human"))

# ST numbers for each human isolate:
HumanST <- ST[(Country==Ctr)&(z)&(Group=="Human")]

positionST <- 1:length(STu)
IST <- matrix(0,Nisolates,length(STu))
for(i in 1:Nisolates){
IST[i,positionST[STu==HumanST[i]]] <- 1 
}

ns<-maxns  # one "unknown source" added

# Overall type frequencies from data (if needed):
OST <- numeric()
for(j in 1:length(STu)){OST[j] <- sum(sourcesST[,j])/sum(sourcesST[]) }
NsourcesST <- sum(sourcesST[])

MCMC <- 10000  # number of MCMC iterations

# VARIABLE STRUCTURES FOR THE MCMC SAMPLER:
numST <- numeric()
ind <- numeric() 
# numbering of ST types:
for(j in 1:length(STu)){  numST[j] <- j }
# build index for ST types of the ith isolate:
prodIST<-numeric()
for(i in 1:Nisolates){
for(j in 1:length(STu)){prodIST[j] <- IST[i,j]*numST[j]}
ind[i] <- sum(prodIST[])
}
prodST <- structure(.Data=c(rep(0,Nisolates*ns*length(STu))),.Dim=c(Nisolates,ns,length(STu)))
STpar0 <- matrix(0,ns,length(STu))
qST<- structure(.Data=c(rep(0,MCMC*ns*length(STu))),.Dim=c(MCMC,ns,length(STu)))
phii <- matrix(0,Nisolates,ns)
phii0 <- matrix(0,Nisolates,ns)
Z <- numeric()
S <- matrix(0,Nisolates,ns) # structure(.Data=c(rep(0,MCMC*Nisolates*ns)),.Dim=c(MCMC,Nisolates,ns))
phi0 <- numeric()
phi <- matrix(0,MCMC,ns) 

###################
# INITIAL VALUES FOR MCMC:
for(i in 1:ns){
for(j in 1:length(STu)){qST[1,i,j] <- 1/length(STu)}
}
phi[1,1:ns] <- rep(1/ns,ns)
Z[1:Nisolates]<-round(runif(Nisolates,0.5/ns,1+0.5/ns)*ns)
for(s in 1:Nisolates){
  for(k in 1:ns){
  S[s,k] <- Z[s]==k   # membership of sth isolate to group k  
  }
}


##############################
##############################
# START MCMC 
# GIBBS SAMPLER *adapted* from 
# Pella & Masuda: Bayesian methods for analysis of stock
# mixtures from genetic characters. Fish. Bull. 99: 151-167 (2001).
h<-0  # control parameter (for avoiding very small values of q ~Dirich, maybe not needed)
b<-0  # sum(humansST)/ns # weighting factor for 'baseline center' = (1/ns)*number.of.human.cases
for(mc in 2:MCMC){ 

# FULL CONDITIONAL DISTRIBUTIONS FOR 
# THE RELATIVE ST TYPE FREQUENCIES IN EACH SOURCE POPULATION:
for(i in 1:ns){   # ns = number of source populations
for(j in 1:length(STu)){ 
  STpar0[i,j] <- rgamma(1,sourcesST[i,j]+sum(S[,i]*IST[,j]>0)+b*sum(sourcesST[,j])/NsourcesST+1/length(STu),1)
}
for(j in 1:length(STu)){
  qST[mc,i,j] <- (STpar0[i,j]+h) / (sum(STpar0[i,])+h*length(STu))
}
}  # end of ns

# FULL CONDITIONAL DISTRIBUTION FOR EACH Z: 
# (Z = group indicators for each isolate) 
for(i in 1:Nisolates){
for(k in 1:ns){ 
phii0[i,k] <- exp( log(phi[mc-1,k])+log(qST[mc,k,ind[i]]) )
}
for(k in 1:ns){
phii[i,k] <- exp(log(phii0[i,k])-log(sum(phii0[i,])) ) # normalize
}
nn <- 1:ns
Z[i] <- min(nn[runif(1)<cumsum(phii[i,])])     # "rcat(1,phii[i,])"
 
for(k in 1:ns){
    S[i,k] <- Z[i]==k   # membership of ith isolate to group k  
    }
} # end Nisolates

# FULL CONDITIONAL DISTRIBUTION FOR phi:
# (phi = proportions of source populations)
for(k in 1:ns){
phi0[k] <- rgamma(1,sum(S[,k])+1/ns,1)
}
for(k in 1:ns){
phi[mc,k] <- phi0[k]/sum(phi0[])
}  
}  # end of MCMC
################################################

burnin <- round(MCMC*0.1)  # set the burn-in period length

par(mfrow=c(ns,1),mar=c(2.5,5.5,2.5,2.5))
for(i in 1:ns){plot(phi[burnin:MCMC,i],xlab="",ylab=sourcenames[i],pch=16,cex=0.4,cex.lab=1.7)}
dev.new()

loST<-matrix(0,ns,length(STu))
upST<-matrix(0,ns,length(STu))

# check the fit of estimated freqs and observed freqs:
# Set the number of picture frames as suitable:
par(mfrow=c(2,1))
errorST <- matrix(NA,ns-1,length(STu))
mqST <- matrix(NA,ns-1,length(STu))
fST <-  matrix(NA,ns-1,length(STu))

for(i in 1:(ns-1)){for(j in 1:length(STu)){
mqST[i,j] <- mean(qST[burnin:MCMC,i,j]) 
fST[i,j] <- sourcesST[i,j]/sum(sourcesST[i,])
loup <- quantile(qST[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
loST[i,j] <- loup[1]; upST[i,j] <- loup[2]
} }
lim <- max(c(max(mqST),max(fST)))
plot(c(0,1),c(0,1),type="l",xlim=c(0,lim*1.1),ylim=c(0,lim*1.1),xlab="Estimated freq",ylab="Observed sample freq",main="ST")
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
plot(density(phi[burnin:MCMC,ix[1]]),col=ix[1],xlim=c(0,1),xlab="Source attribution",ylab="",main=Ctr,lwd=2)
for(i in 2:ns){
points(density(phi[burnin:MCMC,ix[i]]),type="l",col=ix[i],lwd=2)
}
names <- c(sourcenames,rep("Unknown",10))
legend(x="topright",names[ix],text.col=ix)
  #dev.off()
dev.new()

errorsum <- numeric() 
for(i in 1:ns-1){
errorsum[i] <- sum(errorST[i,])
}

ES <- sum(errorsum)   # sum of all squared errors

#################################################
#  Optional:
#################################################
# Compare human freqs with freqs in each source
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Human sample freq",ylab="Source sample freq",main="ST")
for(i in 1:(ns-1)){
points(humansST/sum(humansST),sourcesST[i,]/sum(sourcesST[i,]),col=i,pch=16)
}


