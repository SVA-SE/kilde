
sourcenames=setdiff(unique(Group[Country==Ctr]),"Human")
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
# define "unknown source" (ns+1) 
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

# Get allele numbers for each human isolate:
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

# build index for allele types of the ith isolate:
ind <- matrix(0,Nisolates,7)
prodIASP<-numeric();prodIGLN<-numeric();
prodIGLT<-numeric();prodIGLY<-numeric();
prodIPGM<-numeric();prodITKT<-numeric();
prodIUNC<-numeric();
for(i in 1:Nisolates){
for(j in 1:nat[1]){prodIASP[j] <- IASP[i,j]*j}
ind[i,1] <- sum(prodIASP[])
for(j in 1:nat[2]){prodIGLN[j] <- IGLN[i,j]*j}
ind[i,2] <- sum(prodIGLN[])
for(j in 1:nat[3]){prodIGLT[j] <- IGLT[i,j]*j}
ind[i,3] <- sum(prodIGLT[])
for(j in 1:nat[4]){prodIGLY[j] <- IGLY[i,j]*j}
ind[i,4] <- sum(prodIGLY[])
for(j in 1:nat[5]){prodIPGM[j] <- IPGM[i,j]*j}
ind[i,5] <- sum(prodIPGM[])
for(j in 1:nat[6]){prodITKT[j] <- ITKT[i,j]*j}
ind[i,6] <- sum(prodITKT[])
for(j in 1:nat[7]){prodIUNC[j] <- IUNC[i,j]*j}
ind[i,7] <- sum(prodIUNC[])
}

ns<-maxns  # one "unknown source" added, with data sample zero, or the following:

# If UM=1, set the average sample as "data" for the unknown source:
if(UM==1){   
for(j in 1:nat[1]){sourcesASP[ns,j] <- round(mean(sourcesASP[1:ns-1,j]))}
for(j in 1:nat[2]){sourcesGLN[ns,j] <- round(mean(sourcesGLN[1:ns-1,j]))}
for(j in 1:nat[3]){sourcesGLT[ns,j] <- round(mean(sourcesGLT[1:ns-1,j]))}
for(j in 1:nat[4]){sourcesGLY[ns,j] <- round(mean(sourcesGLY[1:ns-1,j]))}
for(j in 1:nat[5]){sourcesPGM[ns,j] <- round(mean(sourcesPGM[1:ns-1,j]))}
for(j in 1:nat[6]){sourcesTKT[ns,j] <- round(mean(sourcesTKT[1:ns-1,j]))}
for(j in 1:nat[7]){sourcesUNC[ns,j] <- round(mean(sourcesUNC[1:ns-1,j]))}
}

if(UM==2){
# This will set the allele counts for the "unknown source sample" 
# such that it represents the counts drawn from STs that were unique to humans.
for(j in 1:length(STu)){
sourcesASP[ns,positionASP[ASPu==STtable[j,2]]]<-sourcesASP[ns,positionASP[ASPu==STtable[j,2]]]+sum((ST==STu[j])&is.element(STu[j],STuHo)&(Group=="Human")&(Country==Ctr)&z) 
}
for(j in 1:length(STu)){
sourcesGLN[ns,positionGLN[GLNu==STtable[j,3]]]<-sourcesGLN[ns,positionGLN[GLNu==STtable[j,3]]]+sum((ST==STu[j])&is.element(STu[j],STuHo)&(Group=="Human")&(Country==Ctr)&z) 
}
for(j in 1:length(STu)){
sourcesGLT[ns,positionGLT[GLTu==STtable[j,4]]]<-sourcesGLT[ns,positionGLT[GLTu==STtable[j,4]]]+sum((ST==STu[j])&is.element(STu[j],STuHo)&(Group=="Human")&(Country==Ctr)&z) 
}
for(j in 1:length(STu)){
sourcesGLY[ns,positionGLY[GLYu==STtable[j,5]]]<-sourcesGLY[ns,positionGLY[GLYu==STtable[j,5]]]+sum((ST==STu[j])&is.element(STu[j],STuHo)&(Group=="Human")&(Country==Ctr)&z) 
}
for(j in 1:length(STu)){
sourcesPGM[ns,positionPGM[PGMu==STtable[j,6]]]<-sourcesPGM[ns,positionPGM[PGMu==STtable[j,6]]]+sum((ST==STu[j])&is.element(STu[j],STuHo)&(Group=="Human")&(Country==Ctr)&z) 
}
for(j in 1:length(STu)){
sourcesTKT[ns,positionTKT[TKTu==STtable[j,7]]]<-sourcesTKT[ns,positionTKT[TKTu==STtable[j,7]]]+sum((ST==STu[j])&is.element(STu[j],STuHo)&(Group=="Human")&(Country==Ctr)&z) 
}
for(j in 1:length(STu)){
sourcesUNC[ns,positionUNC[UNCu==STtable[j,8]]]<-sourcesUNC[ns,positionUNC[UNCu==STtable[j,8]]]+sum((ST==STu[j])&is.element(STu[j],STuHo)&(Group=="Human")&(Country==Ctr)&z) 
}
}





