
DATA <- read.table("NMDD2015data_vers21_modif.txt",header=TRUE)

attach(DATA)
# use ALL allele types
z <- !is.na(ST)

Ctr = "FI"
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

######################
# make a table of STs and corresponding alleles:
STtable<- matrix(NA,length(STu),8)
for(i in 1:length(STu)){
STtable[i,1] <- STu[i]
STtable[i,2] <- unique(ASP[(ST==STu[i])&!is.na(ST)&Country==Ctr])
STtable[i,3] <- unique(GLN[(ST==STu[i])&!is.na(ST)&Country==Ctr])
STtable[i,4] <- unique(GLT[(ST==STu[i])&!is.na(ST)&Country==Ctr])
STtable[i,5] <- unique(GLY[(ST==STu[i])&!is.na(ST)&Country==Ctr])
STtable[i,6] <- unique(PGM[(ST==STu[i])&!is.na(ST)&Country==Ctr])
STtable[i,7] <- unique(TKT[(ST==STu[i])&!is.na(ST)&Country==Ctr])
STtable[i,8] <- unique(UNC[(ST==STu[i])&!is.na(ST)&Country==Ctr])
}
############################

# find how many different a-types there are, and list them:
ASPu <- sort(unique(STtable[,2]))
GLNu <- sort(unique(STtable[,3]))
GLTu <- sort(unique(STtable[,4]))
GLYu <- sort(unique(STtable[,5]))
PGMu <- sort(unique(STtable[,6]))
TKTu <- sort(unique(STtable[,7]))
UNCu <- sort(unique(STtable[,8]))

nat <- numeric()  # number of all different alleles
nat[1] <- length(ASPu) 
nat[2] <- length(GLNu)
nat[3] <- length(GLTu)
nat[4] <- length(GLYu)
nat[5] <- length(PGMu)
nat[6] <- length(TKTu)
nat[7] <- length(UNCu)

# number (counts) of a-types in each source:
# define "unknown source" (ns+1) with zero observed counts.
# OR: define many unknown sources: ns+K for the DPP model?
maxns <- ns+1
sourcesASP <- matrix(0,maxns,length(ASPu)) 
sourcesGLN <- matrix(0,maxns,length(GLNu))
sourcesGLT <- matrix(0,maxns,length(GLTu))
sourcesGLY <- matrix(0,maxns,length(GLYu))
sourcesPGM <- matrix(0,maxns,length(PGMu))
sourcesTKT <- matrix(0,maxns,length(TKTu))
sourcesUNC <- matrix(0,maxns,length(UNCu))

# number (counts) of such a-types in humans & also found in sources:
humansASP <- numeric()
humansGLN <- numeric()
humansGLT <- numeric()
humansGLY <- numeric()
humansPGM <- numeric()
humansTKT <- numeric()
humansUNC <- numeric()

# sample frequencies (counts) of each a-type
for(i in 1:length(ASPu)){
humansASP[i] <- sum((ASP==ASPu[i])&(Group=="Human")&(Country==Ctr)&z)
for(j in 1:ns){
sourcesASP[j,i] <- sum((ASP==ASPu[i])&(Group==sourcenames[j])&(Country==Ctr)&z) }
}
for(i in 1:length(GLNu)){
humansGLN[i] <- sum((GLN==GLNu[i])&(Group=="Human")&(Country==Ctr)&z)
for(j in 1:ns){
sourcesGLN[j,i] <- sum((GLN==GLNu[i])&(Group==sourcenames[j])&(Country==Ctr)&z)   }
}
for(i in 1:length(GLTu)){
humansGLT[i] <- sum((GLT==GLTu[i])&(Group=="Human")&(Country==Ctr)&z)
for(j in 1:ns){
sourcesGLT[j,i] <- sum((GLT==GLTu[i])&(Group==sourcenames[j])&(Country==Ctr)&z)  }
}
for(i in 1:length(GLYu)){
humansGLY[i] <- sum((GLY==GLYu[i])&(Group=="Human")&(Country==Ctr)&z)
for(j in 1:ns){
sourcesGLY[j,i] <- sum((GLY==GLYu[i])&(Group==sourcenames[j])&(Country==Ctr)&z) }
}
for(i in 1:length(PGMu)){
humansPGM[i] <-sum((PGM==PGMu[i])&(Group=="Human")&(Country==Ctr)&z)
for(j in 1:ns){ 
sourcesPGM[j,i] <- sum((PGM==PGMu[i])&(Group==sourcenames[j])&(Country==Ctr)&z)  }
}
for(i in 1:length(TKTu)){
humansTKT[i] <- sum((TKT==TKTu[i])&(Group=="Human")&(Country==Ctr)&z)
for(j in 1:ns){
sourcesTKT[j,i] <- sum((TKT==TKTu[i])&(Group==sourcenames[j])&(Country==Ctr)&z)   }
}
for(i in 1:length(UNCu)){
humansUNC[i] <- sum((UNC==UNCu[i])&(Group=="Human")&(Country==Ctr)&z)
for(j in 1:ns){
sourcesUNC[j,i] <- sum((UNC==UNCu[i])&(Group==sourcenames[j])&(Country==Ctr)&z)  }
}
# Number of all human isolates:
Nisolates <- sum((Country==Ctr)&(z)&(Group=="Human"))

# Allele numbers for each human isolate:
HumanASP <- ASP[(Country==Ctr)&(z)&(Group=="Human")]
HumanGLN <- GLN[(Country==Ctr)&(z)&(Group=="Human")]
HumanGLT <- GLT[(Country==Ctr)&(z)&(Group=="Human")]
HumanGLY <- GLY[(Country==Ctr)&(z)&(Group=="Human")]
HumanPGM <- PGM[(Country==Ctr)&(z)&(Group=="Human")]
HumanTKT <- TKT[(Country==Ctr)&(z)&(Group=="Human")]
HumanUNC <- UNC[(Country==Ctr)&(z)&(Group=="Human")]

positionASP <- 1:length(ASPu)
IASP <- matrix(0,Nisolates,nat[1])
positionGLN <- 1:length(GLNu)
IGLN <- matrix(0,Nisolates,nat[2])
positionGLT <- 1:length(GLTu)
IGLT <- matrix(0,Nisolates,nat[3])
positionGLY <- 1:length(GLYu)
IGLY <- matrix(0,Nisolates,nat[4])
positionPGM <- 1:length(PGMu)
IPGM <- matrix(0,Nisolates,nat[5])
positionTKT <- 1:length(TKTu)
ITKT <- matrix(0,Nisolates,nat[6])
positionUNC <- 1:length(UNCu)
IUNC <- matrix(0,Nisolates,nat[7])
for(i in 1:Nisolates){
IASP[i,positionASP[ASPu==HumanASP[i]]] <- 1 
IGLN[i,positionGLN[GLNu==HumanGLN[i]]] <- 1 
IGLT[i,positionGLT[GLTu==HumanGLT[i]]] <- 1 
IGLY[i,positionGLY[GLYu==HumanGLY[i]]] <- 1 
IPGM[i,positionPGM[PGMu==HumanPGM[i]]] <- 1 
ITKT[i,positionTKT[TKTu==HumanTKT[i]]] <- 1 
IUNC[i,positionUNC[UNCu==HumanUNC[i]]] <- 1 
}

ns<-maxns  # one "unknown source" added
MCMC <- 600  # number of MCMC iterations

# VARIABLE STRUCTURES FOR THE MCMC SAMPLER:
numASP <- numeric()
numGLN <- numeric()
numGLT <- numeric()
numGLY <- numeric()
numPGM <- numeric()
numTKT <- numeric()
numUNC <- numeric()
ind <- matrix(0,Nisolates,7)
# numbering of allele types:
for(j in 1:nat[1]){  numASP[j] <- j }
for(j in 1:nat[2]){  numGLN[j] <- j }
for(j in 1:nat[3]){  numGLT[j] <- j }
for(j in 1:nat[4]){  numGLY[j] <- j }
for(j in 1:nat[5]){  numPGM[j] <- j }
for(j in 1:nat[6]){  numTKT[j] <- j }
for(j in 1:nat[7]){  numUNC[j] <- j }
# build index for allele types of the ith isolate:
prodIASP<-numeric();prodIGLN<-numeric();
prodIGLT<-numeric();prodIGLY<-numeric();
prodIPGM<-numeric();prodITKT<-numeric();
prodIUNC<-numeric();
for(i in 1:Nisolates){
for(j in 1:nat[1]){prodIASP[j] <- IASP[i,j]*numASP[j]}
ind[i,1] <- sum(prodIASP[])
for(j in 1:nat[2]){prodIGLN[j] <- IGLN[i,j]*numGLN[j]}
ind[i,2] <- sum(prodIGLN[])
for(j in 1:nat[3]){prodIGLT[j] <- IGLT[i,j]*numGLT[j]}
ind[i,3] <- sum(prodIGLT[])
for(j in 1:nat[4]){prodIGLY[j] <- IGLY[i,j]*numGLY[j]}
ind[i,4] <- sum(prodIGLY[])
for(j in 1:nat[5]){prodIPGM[j] <- IPGM[i,j]*numPGM[j]}
ind[i,5] <- sum(prodIPGM[])
for(j in 1:nat[6]){prodITKT[j] <- ITKT[i,j]*numTKT[j]}
ind[i,6] <- sum(prodITKT[])
for(j in 1:nat[7]){prodIUNC[j] <- IUNC[i,j]*numUNC[j]}
ind[i,7] <- sum(prodIUNC[])
}
prodASP <- structure(.Data=c(rep(0,Nisolates*ns*nat[1])),.Dim=c(Nisolates,ns,nat[1]))
prodGLN <- structure(.Data=c(rep(0,Nisolates*ns*nat[2])),.Dim=c(Nisolates,ns,nat[2]))
prodGLT <- structure(.Data=c(rep(0,Nisolates*ns*nat[3])),.Dim=c(Nisolates,ns,nat[3]))
prodGLY <- structure(.Data=c(rep(0,Nisolates*ns*nat[4])),.Dim=c(Nisolates,ns,nat[4]))
prodPGM <- structure(.Data=c(rep(0,Nisolates*ns*nat[5])),.Dim=c(Nisolates,ns,nat[5]))
prodTKT <- structure(.Data=c(rep(0,Nisolates*ns*nat[6])),.Dim=c(Nisolates,ns,nat[6]))
prodUNC <- structure(.Data=c(rep(0,Nisolates*ns*nat[7])),.Dim=c(Nisolates,ns,nat[7]))
ASPpar0 <- matrix(0,ns,nat[1])
GLNpar0 <- matrix(0,ns,nat[2])
GLTpar0 <- matrix(0,ns,nat[3])
GLYpar0 <- matrix(0,ns,nat[4])
PGMpar0 <- matrix(0,ns,nat[5])
TKTpar0 <- matrix(0,ns,nat[6])
UNCpar0 <- matrix(0,ns,nat[7])
qASP<- structure(.Data=c(rep(0,MCMC*ns*nat[1])),.Dim=c(MCMC,ns,nat[1]))
qGLN<- structure(.Data=c(rep(0,MCMC*ns*nat[2])),.Dim=c(MCMC,ns,nat[2]))
qGLT<- structure(.Data=c(rep(0,MCMC*ns*nat[3])),.Dim=c(MCMC,ns,nat[3]))
qGLY<- structure(.Data=c(rep(0,MCMC*ns*nat[4])),.Dim=c(MCMC,ns,nat[4]))
qPGM<- structure(.Data=c(rep(0,MCMC*ns*nat[5])),.Dim=c(MCMC,ns,nat[5]))
qTKT<- structure(.Data=c(rep(0,MCMC*ns*nat[6])),.Dim=c(MCMC,ns,nat[6]))
qUNC<- structure(.Data=c(rep(0,MCMC*ns*nat[7])),.Dim=c(MCMC,ns,nat[7]))
phii <- matrix(0,Nisolates,ns)
phii0 <- matrix(0,Nisolates,ns)
Z <- numeric()  # matrix(0,MCMC,Nisolates)
S <- matrix(0,Nisolates,ns) # structure(.Data=c(rep(0,MCMC*Nisolates*ns)),.Dim=c(MCMC,Nisolates,ns))
phi0 <- numeric()
phi <- matrix(0,MCMC,ns) 

###################
# INITIAL VALUES FOR MCMC:
for(i in 1:ns){
for(j in 1:nat[1]){qASP[1,i,j] <- 1/nat[1]}
for(j in 1:nat[2]){qGLN[1,i,j] <- 1/nat[2]}
for(j in 1:nat[3]){qGLT[1,i,j] <- 1/nat[3]}
for(j in 1:nat[4]){qGLY[1,i,j] <- 1/nat[4]}
for(j in 1:nat[5]){qPGM[1,i,j] <- 1/nat[5]}
for(j in 1:nat[6]){qTKT[1,i,j] <- 1/nat[6]}
for(j in 1:nat[7]){qUNC[1,i,j] <- 1/nat[7]}
}
phi[1,1:ns] <- rep(1/ns,ns)
Z[1:Nisolates]<-round(runif(Nisolates,0.5/ns,1+0.5/ns)*ns)
for(s in 1:Nisolates){
  for(k in 1:ns){
  S[s,k] <- Z[s]==k   # membership of sth isolate to group k  
  }
}
###################
# Overall type frequencies from data (if needed):
OASP <- numeric();OGLN <- numeric();
OGLT <- numeric();OGLY <- numeric();
OPGM <- numeric();OTKT <- numeric();
OUNC <- numeric();
for(j in 1:nat[1]){OASP[j] <- sum(sourcesASP[,j])/sum(sourcesASP[]) }
for(j in 1:nat[2]){OGLN[j] <- sum(sourcesGLN[,j])/sum(sourcesGLN[]) }
for(j in 1:nat[3]){OGLT[j] <- sum(sourcesGLT[,j])/sum(sourcesGLT[]) }
for(j in 1:nat[4]){OGLY[j] <- sum(sourcesGLY[,j])/sum(sourcesGLY[]) }
for(j in 1:nat[5]){OPGM[j] <- sum(sourcesPGM[,j])/sum(sourcesPGM[]) }
for(j in 1:nat[6]){OTKT[j] <- sum(sourcesTKT[,j])/sum(sourcesTKT[]) }
for(j in 1:nat[7]){OUNC[j] <- sum(sourcesUNC[,j])/sum(sourcesUNC[]) }
NsourcesASP <- sum(sourcesASP[]);NsourcesGLN <- sum(sourcesGLN[]);
NsourcesGLT <- sum(sourcesGLT[]);NsourcesGLY <- sum(sourcesGLY[]);
NsourcesPGM <- sum(sourcesPGM[]);NsourcesTKT <- sum(sourcesTKT[]);
NsourcesUNC <- sum(sourcesUNC[]);

##############################
##############################
# START MCMC 
# GIBBS SAMPLER *adapted* from 
# Pella & Masuda: Bayesian methods for analysis of stock
# mixtures from genetic characters. Fish. Bull. 99: 151-167 (2001).
h<-0  # control parameter (for avoiding very small values of q ~Dirich)
for(mc in 2:MCMC){ 

# FULL CONDITIONAL DISTRIBUTIONS FOR 
# THE RELATIVE ALLELE TYPE FREQUENCIES IN EACH SOURCE POPULATION:
# (qASP,qGLN,qGLT,qGLY,qPGM,qTKT,qUNC: these are defined for each loci)
for(i in 1:ns){   # ns = number of source populations
for(j in 1:nat[1]){ 
  for(s in 1:Nisolates){prodASP[s,i,j] <- S[s,i]*IASP[s,j]}
  ASPpar0[i,j] <- rgamma(1,sourcesASP[i,j]+sum(prodASP[,i,j])+sum(sourcesASP[,j])+1/nat[1],1)
  qASP[mc,i,j] <-  (ASPpar0[i,j]+h) / (sum(ASPpar0[i,])+h*nat[1])
}
for(j in 1:nat[2]){ 
  for(s in 1:Nisolates){prodGLN[s,i,j] <- S[s,i]*IGLN[s,j]}
  GLNpar0[i,j] <- rgamma(1,sourcesGLN[i,j]+sum(prodGLN[,i,j])+sum(sourcesGLN[,j])+1/nat[2],1)
  qGLN[mc,i,j] <-  (GLNpar0[i,j]+h) / (sum(GLNpar0[i,])+h*nat[2])
}
for(j in 1:nat[3]){ 
  for(s in 1:Nisolates){prodGLT[s,i,j] <- S[s,i]*IGLT[s,j]}
  GLTpar0[i,j] <- rgamma(1,sourcesGLT[i,j]+sum(prodGLT[,i,j])+sum(sourcesGLT[,j])+1/nat[3],1)
  qGLT[mc,i,j] <-  (GLTpar0[i,j]+h) / (sum(GLTpar0[i,])+h*nat[3])
}
for(j in 1:nat[4]){ 
  for(s in 1:Nisolates){prodGLY[s,i,j] <- S[s,i]*IGLY[s,j]} 
  GLYpar0[i,j] <- rgamma(1,sourcesGLY[i,j]+sum(prodGLY[,i,j])+sum(sourcesGLY[,j])+1/nat[4],1)
  qGLY[mc,i,j] <-  (GLYpar0[i,j]+h) / (sum(GLYpar0[i,])+h*nat[4])
}
for(j in 1:nat[5]){ 
  for(s in 1:Nisolates){prodPGM[s,i,j] <- S[s,i]*IPGM[s,j]}
  PGMpar0[i,j] <- rgamma(1,sourcesPGM[i,j]+sum(prodPGM[,i,j])+sum(sourcesPGM[,j])+1/nat[5],1)
  qPGM[mc,i,j] <-  (PGMpar0[i,j]+h) / (sum(PGMpar0[i,])+h*nat[5])
}
for(j in 1:nat[6]){ 
  for(s in 1:Nisolates){prodTKT[s,i,j] <- S[s,i]*ITKT[s,j]}   
  TKTpar0[i,j] <- rgamma(1,sourcesTKT[i,j]+sum(prodTKT[,i,j])+sum(sourcesTKT[,j])+1/nat[6],1)
  qTKT[mc,i,j] <-  (TKTpar0[i,j]+h) / (sum(TKTpar0[i,])+h*nat[6])
}
for(j in 1:nat[7]){ 
  for(s in 1:Nisolates){prodUNC[s,i,j] <- S[s,i]*IUNC[s,j]}   
  UNCpar0[i,j] <- rgamma(1,sourcesUNC[i,j]+sum(prodUNC[,i,j])+sum(sourcesUNC[,j])+1/nat[7],1)
  qUNC[mc,i,j] <-  (UNCpar0[i,j]+h) / (sum(UNCpar0[i,])+h*nat[7])
}
}  # end of ns

# FULL CONDITIONAL DISTRIBUTION FOR EACH Z: 
# (Z = group indicators for each isolate) 
for(i in 1:Nisolates){
for(k in 1:ns){ 
phii0[i,k] <- exp( log(phi[mc-1,k])+log(qASP[mc,k,ind[i,1]])+log(qGLN[mc,k,ind[i,2]])+log(qGLT[mc,k,ind[i,3]])+log(qGLY[mc,k,ind[i,4]])+log(qPGM[mc,k,ind[i,5]])+log(qTKT[mc,k,ind[i,6]])+log(qUNC[mc,k,ind[i,7]]) )
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
#################

burnin <- 11
loASP<-matrix(0,ns,nat[1]); loGLN<-matrix(0,ns,nat[2])
loGLT<-matrix(0,ns,nat[3]); loGLY<-matrix(0,ns,nat[4])
loPGM<-matrix(0,ns,nat[5]); loTKT<-matrix(0,ns,nat[6])
loUNC<-matrix(0,ns,nat[7])
upASP<-matrix(0,ns,nat[1]); upGLN<-matrix(0,ns,nat[2])
upGLT<-matrix(0,ns,nat[3]); upGLY<-matrix(0,ns,nat[4])
upPGM<-matrix(0,ns,nat[5]); upTKT<-matrix(0,ns,nat[6])
upUNC<-matrix(0,ns,nat[7])

# check the fit of estimated freqs and observed freqs:
# Set the number of picture frames as suitable:
par(mfrow=c(2,4))

mqASP <- matrix(NA,ns-1,nat[1])
for(i in 1:(ns-1)){for(j in 1:nat[1]){
mqASP[i,j] <- mean(qASP[burnin:MCMC,i,j]) 
loup <- quantile(qASP[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
loASP[i,j] <- loup[1]; upASP[i,j] <- loup[2]
} }
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="ASP")
for(i in 1:(ns-1)){
points(mqASP[i,],sourcesASP[i,]/sum(sourcesASP[i,]),col=i,pch=16)
for(j in 1:nat[1]){points(c(loASP[i,j],upASP[i,j]),rep(sourcesASP[i,j]/sum(sourcesASP[i,]),2),col=i,type="l")}
}
mqGLN <- matrix(NA,ns-1,nat[2])
for(i in 1:(ns-1)){for(j in 1:nat[2]){
mqGLN[i,j] <- mean(qGLN[burnin:MCMC,i,j]) 
loup <- quantile(qGLN[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
loGLN[i,j] <- loup[1]; upGLN[i,j] <- loup[2]
} }
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="GLN")
for(i in 1:(ns-1)){
points(mqGLN[i,],sourcesGLN[i,]/sum(sourcesGLN[i,]),col=i,pch=16)
for(j in 1:nat[2]){points(c(loGLN[i,j],upGLN[i,j]),rep(sourcesGLN[i,j]/sum(sourcesGLN[i,]),2),col=i,type="l")}
}
mqGLT <- matrix(NA,ns-1,nat[3])
for(i in 1:(ns-1)){for(j in 1:nat[3]){
mqGLT[i,j] <- mean(qGLT[burnin:MCMC,i,j]) 
loup <- quantile(qGLT[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
loGLT[i,j] <- loup[1]; upGLT[i,j] <- loup[2]
} }
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="GLT")
for(i in 1:(ns-1)){
points(mqGLT[i,],sourcesGLT[i,]/sum(sourcesGLT[i,]),col=i,pch=16)
for(j in 1:nat[3]){points(c(loGLT[i,j],upGLT[i,j]),rep(sourcesGLT[i,j]/sum(sourcesGLT[i,]),2),col=i,type="l")}
}
mqGLY <- matrix(NA,ns-1,nat[4])
for(i in 1:(ns-1)){for(j in 1:nat[4]){
mqGLY[i,j] <- mean(qGLY[burnin:MCMC,i,j]) 
loup <- quantile(qGLY[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
loGLY[i,j] <- loup[1]; upGLY[i,j] <- loup[2]
} }
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="GLY")
for(i in 1:(ns-1)){
points(mqGLY[i,],sourcesGLY[i,]/sum(sourcesGLY[i,]),col=i,pch=16)
for(j in 1:nat[4]){points(c(loGLY[i,j],upGLY[i,j]),rep(sourcesGLY[i,j]/sum(sourcesGLY[i,]),2),col=i,type="l")}
}
mqPGM <- matrix(NA,ns-1,nat[5])
for(i in 1:(ns-1)){for(j in 1:nat[5]){
mqPGM[i,j] <- mean(qPGM[burnin:MCMC,i,j]) 
loup <- quantile(qPGM[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
loPGM[i,j] <- loup[1]; upPGM[i,j] <- loup[2]
} }
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="PGM")
for(i in 1:(ns-1)){
points(mqPGM[i,],sourcesPGM[i,]/sum(sourcesPGM[i,]),col=i,pch=16)
for(j in 1:nat[5]){points(c(loPGM[i,j],upPGM[i,j]),rep(sourcesPGM[i,j]/sum(sourcesPGM[i,]),2),col=i,type="l")}
}
mqTKT <- matrix(NA,ns-1,nat[6])
for(i in 1:(ns-1)){for(j in 1:nat[6]){
mqTKT[i,j] <- mean(qTKT[burnin:MCMC,i,j]) 
loup <- quantile(qTKT[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
loTKT[i,j] <- loup[1]; upTKT[i,j] <- loup[2]
} }
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="TKT")
for(i in 1:(ns-1)){
points(mqTKT[i,],sourcesTKT[i,]/sum(sourcesTKT[i,]),col=i,pch=16)
for(j in 1:nat[6]){points(c(loTKT[i,j],upTKT[i,j]),rep(sourcesTKT[i,j]/sum(sourcesTKT[i,]),2),col=i,type="l")}
}
mqUNC <- matrix(NA,ns-1,nat[7])
for(i in 1:(ns-1)){for(j in 1:nat[7]){
mqUNC[i,j] <- mean(qUNC[burnin:MCMC,i,j]) 
loup <- quantile(qUNC[burnin:MCMC,i,j],c(0.025,0.975),names=FALSE)
loUNC[i,j] <- loup[1]; upUNC[i,j] <- loup[2]
} }
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Estimated freq",ylab="Observed sample freq",main="UNC")
for(i in 1:(ns-1)){
points(mqUNC[i,],sourcesUNC[i,]/sum(sourcesUNC[i,]),col=i,pch=16)
for(j in 1:nat[7]){points(c(loUNC[i,j],upUNC[i,j]),rep(sourcesUNC[i,j]/sum(sourcesUNC[i,]),2),col=i,type="l")}
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


# Compare human freqs with freqs in each source
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Human sample freq",ylab="Source sample freq",main="ASP")
for(i in 1:(ns-1)){
points(humansASP/sum(humansASP),sourcesASP[i,]/sum(sourcesASP[i,]),col=i,pch=16)
}
# Compare human freqs with freqs in each source
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Human sample freq",ylab="Source sample freq",main="GLN")
for(i in 1:(ns-1)){
points(humansGLN/sum(humansGLN),sourcesGLN[i,]/sum(sourcesGLN[i,]),col=i,pch=16)
}
# Compare human freqs with freqs in each source
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Human sample freq",ylab="Source sample freq",main="GLT")
for(i in 1:(ns-1)){
points(humansGLT/sum(humansGLT),sourcesGLT[i,]/sum(sourcesGLT[i,]),col=i,pch=16)
}
# Compare human freqs with freqs in each source
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Human sample freq",ylab="Source sample freq",main="GLY")
for(i in 1:(ns-1)){
points(humansGLY/sum(humansGLY),sourcesGLY[i,]/sum(sourcesGLY[i,]),col=i,pch=16)
}
# Compare human freqs with freqs in each source
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Human sample freq",ylab="Source sample freq",main="PGM")
for(i in 1:(ns-1)){
points(humansPGM/sum(humansPGM),sourcesPGM[i,]/sum(sourcesPGM[i,]),col=i,pch=16)
}
# Compare human freqs with freqs in each source
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Human sample freq",ylab="Source sample freq",main="TKT")
for(i in 1:(ns-1)){
points(humansTKT/sum(humansTKT),sourcesTKT[i,]/sum(sourcesTKT[i,]),col=i,pch=16)
}
# Compare human freqs with freqs in each source
plot(c(0,1),c(0,1),type="l",xlim=c(0,1),ylim=c(0,1),xlab="Human sample freq",ylab="Source sample freq",main="UNC")
for(i in 1:(ns-1)){
points(humansUNC/sum(humansUNC),sourcesUNC[i,]/sum(sourcesUNC[i,]),col=i,pch=16)
}





