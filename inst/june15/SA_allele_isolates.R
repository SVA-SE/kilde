#########
#### ALLELE-model

# Compute attribution isolate by isolate

setwd("C://Users/Jukkis/Documents/Evira/Kampy/SourceAttribution/Mallit/Lopulliset/ToThomas")
DATA <- read.table("NMDD2015data_vers21_modif.txt",header=TRUE)

attach(DATA)
# use ALL allele types
z <- !is.na(ST)

Ctr = "FI"
sourcenames=c("Cattle","Broiler","Turkey","Water") # FI  
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

# number (counts) of a-types:
sourcesASP <- matrix(0,ns+1,length(ASPu)) 
sourcesGLN <- matrix(0,ns+1,length(GLNu))
sourcesGLT <- matrix(0,ns+1,length(GLTu))
sourcesGLY <- matrix(0,ns+1,length(GLYu))
sourcesPGM <- matrix(0,ns+1,length(PGMu))
sourcesTKT <- matrix(0,ns+1,length(TKTu))
sourcesUNC <- matrix(0,ns+1,length(UNCu))

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
#sourcesASP[ns+1,i] <- 0
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
 
library("R2OpenBUGS")

ns <- ns+1  # sources + unknown source

beta <- rep(1,ns)  # Dir prior parameters for sources
data <- list("IASP","IGLN","IGLT","IGLY","IPGM","ITKT","IUNC",
"sourcesASP","sourcesGLN","sourcesGLT","sourcesGLY","sourcesPGM","sourcesTKT","sourcesUNC",
"ns","nat","beta","Nisolates")

parameters <- c("qASP","phi","P","Z")

inits <- function(){list(g0=rep(1,ns))}

res <- bugs(data,inits,parameters,"SA_allele_isolate.txt",n.chains=1,n.burnin=100,n.iter=5100)
attach.bugs(res)

s <- numeric()
for(i in 1:ns){ s[i] <- sd(phi[,i]) }
iix <- sort(s,index.return=TRUE)
attach(iix)

plot(density(phi[,ix[1]],col=ix[1]),xlim=c(0,1),xlab="Source attribution",ylab="",main=Ctr,lwd=2)
for(i in 2:ns){
points(density(phi[,ix[i]]),type="l",col=ix[i],lwd=2)
}
names <- c(sourcenames,"Unknown")
legend(x="topright",names[ix],text.col=ix)

####################################################




