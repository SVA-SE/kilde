#########
#### ALLELE-model


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

# find how many different a-types there are (in all source groups), and list them:
ASPu <- sort(unique(ASP[(Country==Ctr)&(z)&(Group!="Human")]))
GLNu <- sort(unique(GLN[(Country==Ctr)&(z)&(Group!="Human")]))
GLTu <- sort(unique(GLT[(Country==Ctr)&(z)&(Group!="Human")]))
GLYu <- sort(unique(GLY[(Country==Ctr)&(z)&(Group!="Human")]))
PGMu <- sort(unique(PGM[(Country==Ctr)&(z)&(Group!="Human")]))
TKTu <- sort(unique(TKT[(Country==Ctr)&(z)&(Group!="Human")]))
UNCu <- sort(unique(UNC[(Country==Ctr)&(z)&(Group!="Human")]))

# find how many different a-types there are (only in humans), and list them:
ASPuHo <- setdiff(sort(unique(ASP[(Country==Ctr)&(z)&(Group=="Human")])),ASPu)
GLNuHo <- setdiff(sort(unique(GLN[(Country==Ctr)&(z)&(Group=="Human")])),GLNu)
GLTuHo <- setdiff(sort(unique(GLT[(Country==Ctr)&(z)&(Group=="Human")])),GLTu)
GLYuHo <- setdiff(sort(unique(GLY[(Country==Ctr)&(z)&(Group=="Human")])),GLYu)
PGMuHo <- setdiff(sort(unique(PGM[(Country==Ctr)&(z)&(Group=="Human")])),PGMu)
TKTuHo <- setdiff(sort(unique(TKT[(Country==Ctr)&(z)&(Group=="Human")])),TKTu)
UNCuHo <- setdiff(sort(unique(UNC[(Country==Ctr)&(z)&(Group=="Human")])),UNCu)


nat <- numeric()  # number of all different alleles (in all sources)
nat[1] <- length(ASPu) 
nat[2] <- length(GLNu)
nat[3] <- length(GLTu)
nat[4] <- length(GLYu)
nat[5] <- length(PGMu)
nat[6] <- length(TKTu)
nat[7] <- length(UNCu)

natHo <- numeric()  # number of all different alleles (only found in humans)
natHo[1] <- length(ASPuHo) 
natHo[2] <- length(GLNuHo)
natHo[3] <- length(GLTuHo)
natHo[4] <- length(GLYuHo)
natHo[5] <- length(PGMuHo)
natHo[6] <- length(TKTuHo)
natHo[7] <- length(UNCuHo)

# number (counts) of a-types found in sources:
sourcesASP <- matrix(NA,ns,length(ASPu)) 
sourcesGLN <- matrix(NA,ns,length(GLNu))
sourcesGLT <- matrix(NA,ns,length(GLTu))
sourcesGLY <- matrix(NA,ns,length(GLYu))
sourcesPGM <- matrix(NA,ns,length(PGMu))
sourcesTKT <- matrix(NA,ns,length(TKTu))
sourcesUNC <- matrix(NA,ns,length(UNCu))

# number (counts) of such a-types in humans & also found in sources:
humansASP <- numeric()
humansGLN <- numeric()
humansGLT <- numeric()
humansGLY <- numeric()
humansPGM <- numeric()
humansTKT <- numeric()
humansUNC <- numeric()

# sample frequencies (counts) of each a-type (the a-types found in sources)
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

# number (counts) of a-types found only in humans:
humansASPo <- numeric()
humansGLNo <- numeric()
humansGLTo <- numeric()
humansGLYo <- numeric()
humansPGMo <- numeric()
humansTKTo <- numeric()
humansUNCo <- numeric()

# sample frequencies (counts) of each a-type (the a-types found only in humans):
for(i in 1:length(ASPuHo)){
humansASPo[i] <- sum((ASP==ASPuHo[i])&(Group=="Human")&(Country==Ctr)&z)
}
for(i in 1:length(GLNuHo)){
humansGLNo[i] <- sum((GLN==GLNuHo[i])&(Group=="Human")&(Country==Ctr)&z)
}
for(i in 1:length(GLTuHo)){
humansGLTo[i] <- sum((GLT==GLTuHo[i])&(Group=="Human")&(Country==Ctr)&z)
}
for(i in 1:length(GLYuHo)){
humansGLYo[i] <- sum((GLY==GLYuHo[i])&(Group=="Human")&(Country==Ctr)&z)
}
for(i in 1:length(PGMuHo)){
humansPGMo[i] <- sum((PGM==PGMuHo[i])&(Group=="Human")&(Country==Ctr)&z)
}
for(i in 1:length(TKTuHo)){
humansTKTo[i] <- sum((TKT==TKTuHo[i])&(Group=="Human")&(Country==Ctr)&z)
}
for(i in 1:length(UNCuHo)){
humansUNCo[i] <- sum((UNC==UNCuHo[i])&(Group=="Human")&(Country==Ctr)&z)
}


# number of all human isolates:
Nhumans <- sum((Country==Ctr)&(z)&(Group=="Human"))

# number of human isolates with an a-type not found in sources:
# (out of a total number of Nhumans). 
# To be used in BUGS for estimating proportion of "unknown sources"
NASPHo <- sum(humansASPo)
NGLNHo <- sum(humansGLNo)
NGLTHo <- sum(humansGLTo)
NGLYHo <- sum(humansGLYo)
NPGMHo <- sum(humansPGMo)
NTKTHo <- sum(humansTKTo)
NUNCHo <- sum(humansUNCo)


# "N***" is to be used as sample size in the binomial model for 
# "mutations" ("novel" a-type or "distinct" a-type: not found in other sources) 
# per source (for each locus):
NASP <- numeric();NGLN <- numeric();
NGLT <- numeric();NGLY <- numeric();
NPGM <- numeric();NTKT <- numeric();NUNC <- numeric();

for(i in 1:ns){
NASP[i] <- sum( (Country==Ctr)&(Group==sourcenames[i])&z  ) 
NGLN[i] <- sum( (Country==Ctr)&(Group==sourcenames[i])&z  ) 
NGLT[i] <- sum( (Country==Ctr)&(Group==sourcenames[i])&z  ) 
NGLY[i] <- sum( (Country==Ctr)&(Group==sourcenames[i])&z  ) 
NPGM[i] <- sum( (Country==Ctr)&(Group==sourcenames[i])&z  ) 
NTKT[i] <- sum( (Country==Ctr)&(Group==sourcenames[i])&z  ) 
NUNC[i] <- sum( (Country==Ctr)&(Group==sourcenames[i])&z  ) 
}

# BELOW: loop over sources, list allele types in each source

maxlength <- numeric(); ASPle <-numeric(); GLNle <-numeric();
GLTle <- numeric(); GLYle <- numeric(); PGMle <- numeric();
TKTle <- numeric(); UNCle <- numeric()
for(i in 1:ns){ 
 ASPle[i]<-length(sort(unique(ASP[(Group==sourcenames[i])&(Country==Ctr)&z])))
 GLNle[i]<-length(sort(unique(GLN[(Group==sourcenames[i])&(Country==Ctr)&z])))
 GLTle[i]<-length(sort(unique(GLT[(Group==sourcenames[i])&(Country==Ctr)&z])))
 GLYle[i]<-length(sort(unique(GLY[(Group==sourcenames[i])&(Country==Ctr)&z])))
 PGMle[i]<-length(sort(unique(PGM[(Group==sourcenames[i])&(Country==Ctr)&z])))
 TKTle[i]<-length(sort(unique(TKT[(Group==sourcenames[i])&(Country==Ctr)&z])))
 UNCle[i]<-length(sort(unique(UNC[(Group==sourcenames[i])&(Country==Ctr)&z])))
}
maxlength[1] <- max(ASPle);maxlength[2] <- max(GLNle);
maxlength[3] <- max(GLTle);maxlength[4] <- max(GLYle);
maxlength[5] <- max(PGMle);maxlength[6] <- max(TKTle);maxlength[7] <- max(UNCle);  

ASPusource <- matrix(NA,ns,maxlength[1])
GLNusource <- matrix(NA,ns,maxlength[2])
GLTusource <- matrix(NA,ns,maxlength[3])
GLYusource <- matrix(NA,ns,maxlength[4])
PGMusource <- matrix(NA,ns,maxlength[5])
TKTusource <- matrix(NA,ns,maxlength[6])
UNCusource <- matrix(NA,ns,maxlength[7])

# the list of found a-types for each source, each locus:
# ***usource[] contains listing of allele types found in each source
for(i in 1:ns){  
ASPusource[i,1:ASPle[i]] <-sort(unique(ASP[(Group==sourcenames[i])&(Country==Ctr)&z])) 
GLNusource[i,1:GLNle[i]] <-sort(unique(GLN[(Group==sourcenames[i])&(Country==Ctr)&z])) 
GLTusource[i,1:GLTle[i]] <-sort(unique(GLT[(Group==sourcenames[i])&(Country==Ctr)&z])) 
GLYusource[i,1:GLYle[i]] <-sort(unique(GLY[(Group==sourcenames[i])&(Country==Ctr)&z])) 
PGMusource[i,1:PGMle[i]] <-sort(unique(PGM[(Group==sourcenames[i])&(Country==Ctr)&z])) 
TKTusource[i,1:TKTle[i]] <-sort(unique(TKT[(Group==sourcenames[i])&(Country==Ctr)&z])) 
UNCusource[i,1:UNCle[i]] <-sort(unique(UNC[(Group==sourcenames[i])&(Country==Ctr)&z])) 
}

ASPnovel <- matrix(NA,ns,nat[1]);GLNnovel <- matrix(NA,ns,nat[2]);GLTnovel <- matrix(NA,ns,nat[3]);
GLYnovel <- matrix(NA,ns,nat[4]);PGMnovel <- matrix(NA,ns,nat[5]);TKTnovel <- matrix(NA,ns,nat[6]);
UNCnovel <- matrix(NA,ns,nat[7]); 
aspnn <- numeric(); glnnn <- numeric(); gltnn <- numeric(); 
glynn <- numeric(); pgmnn <- numeric(); tktnn <- numeric(); uncnn <- numeric(); 
  
# Find out list of 'novel types' in a source (that were not found in other sources)
for(i in 1:ns){ 
aspnn[i]<-length(setdiff(ASPusource[i,1:ASPle[i]],sort(unique(ASP[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z]))))
if(aspnn[i]>0){
ASPnovel[i,1:aspnn[i]] <- setdiff(ASPusource[i,1:ASPle[i]],sort(unique(ASP[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z])))
}
glnnn[i]<-length(setdiff(GLNusource[i,1:GLNle[i]],sort(unique(GLN[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z]))))
if(glnnn[i]>0){
GLNnovel[i,1:glnnn[i]] <- setdiff(GLNusource[i,1:GLNle[i]],sort(unique(GLN[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z])))
}
gltnn[i]<-length(setdiff(GLTusource[i,1:GLTle[i]],sort(unique(GLT[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z]))))
if(gltnn[i]>0){
GLTnovel[i,1:gltnn[i]] <- setdiff(GLTusource[i,1:GLTle[i]],sort(unique(GLT[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z])))
}
glynn[i]<-length(setdiff(GLYusource[i,1:GLYle[i]],sort(unique(GLY[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z]))))
if(glynn[i]>0){
GLYnovel[i,1:glynn[i]] <- setdiff(GLYusource[i,1:GLYle[i]],sort(unique(GLY[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z])))
}
pgmnn[i]<-length(setdiff(PGMusource[i,1:PGMle[i]],sort(unique(PGM[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z]))))
if(pgmnn[i]>0){
PGMnovel[i,1:pgmnn[i]] <- setdiff(PGMusource[i,1:PGMle[i]],sort(unique(PGM[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z])))
}
tktnn[i]<-length(setdiff(TKTusource[i,1:TKTle[i]],sort(unique(TKT[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z]))))
if(tktnn[i]>0){
TKTnovel[i,1:tktnn[i]] <- setdiff(TKTusource[i,1:TKTle[i]],sort(unique(TKT[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z])))
}
uncnn[i]<-length(setdiff(UNCusource[i,1:UNCle[i]],sort(unique(UNC[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z]))))
if(uncnn[i]>0){
UNCnovel[i,1:uncnn[i]] <- setdiff(UNCusource[i,1:UNCle[i]],sort(unique(UNC[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z])))
}
}

# count the frequencies (numbers) of the novel types:
sourcesASPnovel=matrix(0,ns,nat[1])
sourcesGLNnovel=matrix(0,ns,nat[2])
sourcesGLTnovel=matrix(0,ns,nat[3])
sourcesGLYnovel=matrix(0,ns,nat[4])
sourcesPGMnovel=matrix(0,ns,nat[5])
sourcesTKTnovel=matrix(0,ns,nat[6])
sourcesUNCnovel=matrix(0,ns,nat[7])

for(j in 1:ns){
 for(i in 1:nat[1]){
 sourcesASPnovel[j,i] <- sourcesASP[j,i]*is.element(ASPu[i],ASPnovel[j,1:aspnn[j]]) }
 for(i in 1:nat[2]){
 sourcesGLNnovel[j,i] <- sourcesGLN[j,i]*is.element(GLNu[i],GLNnovel[j,1:glnnn[j]]) }
 for(i in 1:nat[3]){
 sourcesGLTnovel[j,i] <- sourcesGLT[j,i]*is.element(GLTu[i],GLTnovel[j,1:gltnn[j]]) }
 for(i in 1:nat[4]){
 sourcesGLYnovel[j,i] <- sourcesGLY[j,i]*is.element(GLYu[i],GLYnovel[j,1:glynn[j]]) }
 for(i in 1:nat[5]){
 sourcesPGMnovel[j,i] <- sourcesPGM[j,i]*is.element(PGMu[i],PGMnovel[j,1:pgmnn[j]]) }
 for(i in 1:nat[6]){
 sourcesTKTnovel[j,i] <- sourcesTKT[j,i]*is.element(TKTu[i],TKTnovel[j,1:tktnn[j]]) }
 for(i in 1:nat[7]){
 sourcesUNCnovel[j,i] <- sourcesUNC[j,i]*is.element(UNCu[i],UNCnovel[j,1:uncnn[j]]) }
 }   

# number of all 'novel types' per source for BUGS-model:  
# (a-types that were not in other sources)
NASPnovel <- numeric()
for(i in 1:ns){ NASPnovel[i] <- sum(sourcesASPnovel[i,]) }
NGLNnovel <- numeric()
for(i in 1:ns){ NGLNnovel[i] <- sum(sourcesGLNnovel[i,]) }
NGLTnovel <- numeric()
for(i in 1:ns){ NGLTnovel[i] <- sum(sourcesGLTnovel[i,]) }
NGLYnovel <- numeric()
for(i in 1:ns){ NGLYnovel[i] <- sum(sourcesGLYnovel[i,]) }
NPGMnovel <- numeric()
for(i in 1:ns){ NPGMnovel[i] <- sum(sourcesPGMnovel[i,]) }
NTKTnovel <- numeric()
for(i in 1:ns){ NTKTnovel[i] <- sum(sourcesTKTnovel[i,]) }
NUNCnovel <- numeric()
for(i in 1:ns){ NUNCnovel[i] <- sum(sourcesUNCnovel[i,]) }

### find out what a-types there are in a source that were found at least 
# in some other sources:
ASPnco <- numeric(); GLNnco <- numeric(); GLTnco <- numeric(); GLYnco <- numeric(); 
PGMnco <- numeric(); TKTnco <- numeric(); UNCnco <- numeric();    
for(i in 1:ns){
ASPnco[i] <- length(setdiff(ASPusource[i,],ASPnovel[i,]))
GLNnco[i] <- length(setdiff(GLNusource[i,],GLNnovel[i,]))
GLTnco[i] <- length(setdiff(GLTusource[i,],GLTnovel[i,]))
GLYnco[i] <- length(setdiff(GLYusource[i,],GLYnovel[i,]))
PGMnco[i] <- length(setdiff(PGMusource[i,],PGMnovel[i,]))
TKTnco[i] <- length(setdiff(TKTusource[i,],TKTnovel[i,]))
UNCnco[i] <- length(setdiff(UNCusource[i,],UNCnovel[i,]))
}
ASPcommon <- matrix(NA,ns,max(ASPnco)); GLNcommon <- matrix(NA,ns,max(GLNnco));
GLTcommon <- matrix(NA,ns,max(GLTnco)); GLYcommon <- matrix(NA,ns,max(GLYnco));
PGMcommon <- matrix(NA,ns,max(PGMnco)); TKTcommon <- matrix(NA,ns,max(TKTnco));
UNCcommon <- matrix(NA,ns,max(UNCnco)); 
for(i in 1:ns){
ASPcommon[i,1:ASPnco[i]] <- setdiff(ASPusource[i,],ASPnovel[i,])
GLNcommon[i,1:GLNnco[i]] <- setdiff(GLNusource[i,],GLNnovel[i,])
GLTcommon[i,1:GLTnco[i]] <- setdiff(GLTusource[i,],GLTnovel[i,])
GLYcommon[i,1:GLYnco[i]] <- setdiff(GLYusource[i,],GLYnovel[i,])
PGMcommon[i,1:PGMnco[i]] <- setdiff(PGMusource[i,],PGMnovel[i,])
TKTcommon[i,1:TKTnco[i]] <- setdiff(TKTusource[i,],TKTnovel[i,])
UNCcommon[i,1:UNCnco[i]] <- setdiff(UNCusource[i,],UNCnovel[i,])
}
migrateASP <- matrix(NA,ns,nat[1]); migrateGLN <- matrix(NA,ns,nat[2]); 
migrateGLT <- matrix(NA,ns,nat[3]); migrateGLY <- matrix(NA,ns,nat[4]); 
migratePGM <- matrix(NA,ns,nat[5]); migrateTKT <- matrix(NA,ns,nat[6]); 
migrateUNC <- matrix(NA,ns,nat[7]); 

# count frequencies of "migrated" types per source:
for(j in 1:ns){
for(i in 1:nat[1]){
migrateASP[j,i] <- sourcesASP[j,i]*is.element(ASPu[i],ASPcommon[j,])
}
for(i in 1:nat[2]){
migrateGLN[j,i] <- sourcesGLN[j,i]*is.element(GLNu[i],GLNcommon[j,])
}
for(i in 1:nat[3]){
migrateGLT[j,i] <- sourcesGLT[j,i]*is.element(GLTu[i],GLTcommon[j,])
}
for(i in 1:nat[4]){
migrateGLY[j,i] <- sourcesGLY[j,i]*is.element(GLYu[i],GLYcommon[j,])
}
for(i in 1:nat[5]){
migratePGM[j,i] <- sourcesPGM[j,i]*is.element(PGMu[i],PGMcommon[j,])
}
for(i in 1:nat[6]){
migrateTKT[j,i] <- sourcesTKT[j,i]*is.element(TKTu[i],TKTcommon[j,])
}
for(i in 1:nat[7]){
migrateUNC[j,i] <- sourcesUNC[j,i]*is.element(UNCu[i],UNCcommon[j,])
}
}

library("R2OpenBUGS")

# the migration-mutation-unknownsource -model:

beta <- rep(1,ns) # Dir prior parameters for sources
alpha <- structure(.Data=rep(1,ns*ns),.Dim=c(ns,ns)) # Dir prior parameters for migration probabilities
data <- list("humansASP","humansGLN","humansGLT","humansGLY","humansPGM","humansTKT","humansUNC",
"migrateASP","migrateGLN","migrateGLT","migrateGLY","migratePGM","migrateTKT","migrateUNC",
"sourcesASPnovel","sourcesGLNnovel","sourcesGLTnovel","sourcesGLYnovel","sourcesPGMnovel","sourcesTKTnovel","sourcesUNCnovel",
"sourcesASP","sourcesGLN","sourcesGLT","sourcesGLY","sourcesPGM","sourcesTKT","sourcesUNC",
"ns","nat","beta","alpha",
"NASP","NGLN","NGLT","NGLY","NPGM","NTKT","NUNC",
"NASPnovel","NGLNnovel","NGLTnovel","NGLYnovel","NPGMnovel","NTKTnovel","NUNCnovel",
"NASPHo","NGLNHo","NGLTHo","NGLYHo","NPGMHo","NTKTHo","NUNCHo","Nhumans")

parameters <- c("phi","etaASP","khi","pmuta","pother")
inits <- function(){list(g0=rep(1,ns),
pmuta=structure(.Data=rep(0.5,ns*7),.Dim=c(ns,7)),
h0=structure(.Data=rep(1,ns*ns),.Dim=c(ns,ns)))}

resnewFI <- bugs(data,inits,parameters,"SA_allelemodel2New.txt",n.chains=1,n.burnin=100,n.iter=5100)
attach.bugs(resnewFI)


mm <- numeric();for(i in 1:ns){mm[i]<-mean(phi[,i])}
mmo <- order(mm)
ch <- mmo[ns]

par(mfrow=c(ns,1))
for(i in 1:ns){
plot(density(phi[,i]),main="marginal posterior",xlab=sourcenames[i],xlim=c(0,1));points(mean(phi[,i]),0,col="red",pch=16) 
if(i==ch){points(mean(phi[,i]),0,col="blue",pch=16,cex=2);text(mean(phi[,i]),1,"MaxP",cex=1.3)}
}


# For comparison: the simple model with no genetic parameters:

beta <- rep(1,ns)  # Dir prior parameters for sources
data <- list("humansASP","humansGLN","humansGLT","humansGLY","humansPGM","humansTKT","humansUNC",
"sourcesASP","sourcesGLN","sourcesGLT","sourcesGLY","sourcesPGM","sourcesTKT","sourcesUNC",
"ns","nat","beta")

parameters <- c("qASP","phi","etaASP","khi")
inits <- function(){list(g0=rep(1,ns))}

res1 <- bugs(data,inits,parameters,"SA_allelemodel1.txt",n.chains=1,n.burnin=100,n.iter=5100)
attach.bugs(res1)

mm <- numeric();for(i in 1:ns){mm[i]<-mean(phi[,i])}
mmo <- order(mm)
ch <- mmo[ns]

par(mfrow=c(ns,1))
for(i in 1:ns){
plot(density(phi[,i]),main="marginal posterior",xlab=sourcenames[i],xlim=c(0,1));points(mean(phi[,i]),0,col="red",pch=16) 
if(i==ch){points(mean(phi[,i]),0,col="blue",pch=16,cex=2);text(mean(phi[,i]),1,"MaxP",cex=1.3)}
}








