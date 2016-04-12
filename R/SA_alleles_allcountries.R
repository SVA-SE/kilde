#########
#### ALLELE-model

DATA <- read.table("../inst/extdata/NMDD2015data_vers21_modif.txt", header=TRUE)

attach(DATA)
# use ALL allele types
z <- !is.na(ST)

Ctr = "FI"
sourcenames=c("Cattle","Broiler","Turkey","Water")
ns <- length(sourcenames)

ASPu <- sort(unique(ASP[(Country==Ctr)&(z)])) # find how many different a-types there are, and list them.
GLNu <- sort(unique(GLN[(Country==Ctr)&(z)]))
GLTu <- sort(unique(GLT[(Country==Ctr)&(z)]))
GLYu <- sort(unique(GLY[(Country==Ctr)&(z)]))
PGMu <- sort(unique(PGM[(Country==Ctr)&(z)]))
TKTu <- sort(unique(TKT[(Country==Ctr)&(z)]))
UNCu <- sort(unique(UNC[(Country==Ctr)&(z)]))

nat <- numeric()
nat[1] <- length(ASPu)  # number of all different alleles
nat[2] <- length(GLNu)
nat[3] <- length(GLTu)
nat[4] <- length(GLYu)
nat[5] <- length(PGMu)
nat[6] <- length(TKTu)
nat[7] <- length(UNCu)

sourcesASP <- matrix(NA,ns,length(ASPu)) 
sourcesGLN <- matrix(NA,ns,length(GLNu))
sourcesGLT <- matrix(NA,ns,length(GLTu))
sourcesGLY <- matrix(NA,ns,length(GLYu))
sourcesPGM <- matrix(NA,ns,length(PGMu))
sourcesTKT <- matrix(NA,ns,length(TKTu))
sourcesUNC <- matrix(NA,ns,length(UNCu))

humansASP <- numeric()
humansGLN <- numeric()
humansGLT <- numeric()
humansGLY <- numeric()
humansPGM <- numeric()
humansTKT <- numeric()
humansUNC <- numeric()

for(i in 1:length(ASPu)){
humansASP[i] <- sum((ASP[z]==ASPu[i])&(Group[z]=="Human")&(Country[z]==Ctr))
for(j in 1:ns){
sourcesASP[j,i] <- sum((ASP[z]==ASPu[i])&(Group[z]==sourcenames[j])&(Country[z]==Ctr)) }
}
for(i in 1:length(GLNu)){
humansGLN[i] <- sum((GLN[z]==GLNu[i])&(Group[z]=="Human")&(Country[z]==Ctr))
for(j in 1:ns){
sourcesGLN[j,i] <- sum((GLN[z]==GLNu[i])&(Group[z]==sourcenames[j])&(Country[z]==Ctr))   }
}
for(i in 1:length(GLTu)){
humansGLT[i] <- sum((GLT[z]==GLTu[i])&(Group[z]=="Human")&(Country[z]==Ctr))
for(j in 1:ns){
sourcesGLT[j,i] <- sum((GLT[z]==GLTu[i])&(Group[z]==sourcenames[j])&(Country[z]==Ctr))  }
}
for(i in 1:length(GLYu)){
humansGLY[i] <- sum((GLY[z]==GLYu[i])&(Group[z]=="Human")&(Country[z]==Ctr))
for(j in 1:ns){
sourcesGLY[j,i] <- sum((GLY[z]==GLYu[i])&(Group[z]==sourcenames[j])&(Country[z]==Ctr)) }
}
for(i in 1:length(PGMu)){
humansPGM[i] <-sum((PGM[z]==PGMu[i])&(Group[z]=="Human")&(Country[z]==Ctr))
for(j in 1:ns){ 
sourcesPGM[j,i] <- sum((PGM[z]==PGMu[i])&(Group[z]==sourcenames[j])&(Country[z]==Ctr))  }
}
for(i in 1:length(TKTu)){
humansTKT[i] <- sum((TKT[z]==TKTu[i])&(Group[z]=="Human")&(Country[z]==Ctr))
for(j in 1:ns){
sourcesTKT[j,i] <- sum((TKT[z]==TKTu[i])&(Group[z]==sourcenames[j])&(Country[z]==Ctr))   }
}
for(i in 1:length(UNCu)){
humansUNC[i] <- sum((UNC[z]==UNCu[i])&(Group[z]=="Human")&(Country[z]==Ctr))
for(j in 1:ns){
sourcesUNC[j,i] <- sum((UNC[z]==UNCu[i])&(Group[z]==sourcenames[j])&(Country[z]==Ctr))  }
}

# "N***" is to be used as sample size in the binomial model for mutations (novel types) 
# per source (for each locus):
NASP <- numeric();NGLN <- numeric();
NGLT <- numeric();NGLY <- numeric();
NPGM <- numeric();NTKT <- numeric();NUNC <- numeric();

for(i in 1:ns){
NASP[i] <- sum( (Country[z]==Ctr)&(Group[z]==sourcenames[i])  ) 
NGLN[i] <- sum( (Country[z]==Ctr)&(Group[z]==sourcenames[i])  ) 
NGLT[i] <- sum( (Country[z]==Ctr)&(Group[z]==sourcenames[i])  ) 
NGLY[i] <- sum( (Country[z]==Ctr)&(Group[z]==sourcenames[i])  ) 
NPGM[i] <- sum( (Country[z]==Ctr)&(Group[z]==sourcenames[i])  ) 
NTKT[i] <- sum( (Country[z]==Ctr)&(Group[z]==sourcenames[i])  ) 
NUNC[i] <- sum( (Country[z]==Ctr)&(Group[z]==sourcenames[i])  ) 
}

# BELOW: loop over sources, list allele types in each source

maxlength <- numeric(); ASPle <-numeric(); GLNle <-numeric();
GLTle <- numeric(); GLYle <- numeric(); PGMle <- numeric();
TKTle <- numeric(); UNCle <- numeric()
for(i in 1:ns){ 
 ASPle[i]<-length(sort(unique(ASP[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])))
 GLNle[i]<-length(sort(unique(GLN[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])))
 GLTle[i]<-length(sort(unique(GLT[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])))
 GLYle[i]<-length(sort(unique(GLY[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])))
 PGMle[i]<-length(sort(unique(PGM[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])))
 TKTle[i]<-length(sort(unique(TKT[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])))
 UNCle[i]<-length(sort(unique(UNC[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])))
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

# ***usource[] contains listing of allele types found in each source
for(i in 1:ns){
ASPusource[i,1:ASPle[i]] <-sort(unique(ASP[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])) 
GLNusource[i,1:GLNle[i]] <-sort(unique(GLN[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])) 
GLTusource[i,1:GLTle[i]] <-sort(unique(GLT[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])) 
GLYusource[i,1:GLYle[i]] <-sort(unique(GLY[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])) 
PGMusource[i,1:PGMle[i]] <-sort(unique(PGM[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])) 
TKTusource[i,1:TKTle[i]] <-sort(unique(TKT[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])) 
UNCusource[i,1:UNCle[i]] <-sort(unique(UNC[(Group[z]==sourcenames[i])&(Country[z]==Ctr)])) 
}

ASPnovel <- matrix(NA,ns,max(nat));GLNnovel <- matrix(NA,ns,max(nat));GLTnovel <- matrix(NA,ns,max(nat));
GLYnovel <- matrix(NA,ns,max(nat));PGMnovel <- matrix(NA,ns,max(nat));TKTnovel <- matrix(NA,ns,max(nat));
UNCnovel <- matrix(NA,ns,max(nat)); 
aspnn <- numeric(); glnnn <- numeric(); gltnn <- numeric(); 
glynn <- numeric(); pgmnn <- numeric(); tktnn <- numeric(); uncnn <- numeric(); 
  
# Find out 'novel types' in a source (that were not found in other sources)
for(i in 1:ns){ 
aspnn[i]<-length(setdiff(ASPusource[i,1:ASPle[i]],sort(unique(ASP[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)]))))
if(aspnn[i]>0){
ASPnovel[i,1:aspnn[i]] <- setdiff(ASPusource[i,1:ASPle[i]],sort(unique(ASP[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)])))
}
glnnn[i]<-length(setdiff(GLNusource[i,1:GLNle[i]],sort(unique(GLN[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)]))))
if(glnnn[i]>0){
GLNnovel[i,1:glnnn[i]] <- setdiff(GLNusource[i,1:GLNle[i]],sort(unique(GLN[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)])))
}
gltnn[i]<-length(setdiff(GLTusource[i,1:GLTle[i]],sort(unique(GLT[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)]))))
if(gltnn[i]>0){
GLTnovel[i,1:gltnn[i]] <- setdiff(GLTusource[i,1:GLTle[i]],sort(unique(GLT[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)])))
}
glynn[i]<-length(setdiff(GLYusource[i,1:GLYle[i]],sort(unique(GLY[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)]))))
if(glynn[i]>0){
GLYnovel[i,1:glynn[i]] <- setdiff(GLYusource[i,1:GLYle[i]],sort(unique(GLY[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)])))
}
pgmnn[i]<-length(setdiff(PGMusource[i,1:PGMle[i]],sort(unique(PGM[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)]))))
if(pgmnn[i]>0){
PGMnovel[i,1:pgmnn[i]] <- setdiff(PGMusource[i,1:PGMle[i]],sort(unique(PGM[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)])))
}
tktnn[i]<-length(setdiff(TKTusource[i,1:TKTle[i]],sort(unique(TKT[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)]))))
if(tktnn[i]>0){
TKTnovel[i,1:tktnn[i]] <- setdiff(TKTusource[i,1:TKTle[i]],sort(unique(TKT[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)])))
}
uncnn[i]<-length(setdiff(UNCusource[i,1:UNCle[i]],sort(unique(UNC[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)]))))
if(uncnn[i]>0){
UNCnovel[i,1:uncnn[i]] <- setdiff(UNCusource[i,1:UNCle[i]],sort(unique(UNC[((Group[z]!=sourcenames[i])&(Group[z]!="Human"))&(Country[z]==Ctr)])))
}
}
# number of 'novel types' per source for BUGS-model:  (types that were not in other sources)
NASPnovel <- numeric()
for(i in 1:ns){ NASPnovel[i] <- aspnn[i] }
NGLNnovel <- numeric()
for(i in 1:ns){ NGLNnovel[i] <- glnnn[i] }
NGLTnovel <- numeric()
for(i in 1:ns){ NGLTnovel[i] <- gltnn[i] }
NGLYnovel <- numeric()
for(i in 1:ns){ NGLYnovel[i] <- glynn[i] }
NPGMnovel <- numeric()
for(i in 1:ns){ NPGMnovel[i] <- pgmnn[i] }
NTKTnovel <- numeric()
for(i in 1:ns){ NTKTnovel[i] <- tktnn[i] }
NUNCnovel <- numeric()
for(i in 1:ns){ NUNCnovel[i] <- uncnn[i] }


#### Human allele data, find the number of novel types (that are not in sources):
Nhum <- numeric(); Nhumnovel <- numeric()
ASPuhum <- sort(unique(ASP[(Country[z]==Ctr)&(Group[z]=="Human")]))
Nhumnovel[1] <- length(setdiff(ASPuhum,sort(unique(ASP[((Group[z]!="Human"))&(Country[z]==Ctr)]))))
GLNuhum <- sort(unique(GLN[(Country[z]==Ctr)&(Group[z]=="Human")])) 
Nhumnovel[2] <- length(setdiff(GLNuhum,sort(unique(GLN[((Group[z]!="Human"))&(Country[z]==Ctr)]))))
GLTuhum <- sort(unique(GLT[(Country[z]==Ctr)&(Group[z]=="Human")])) 
Nhumnovel[3] <- length(setdiff(GLTuhum,sort(unique(GLT[((Group[z]!="Human"))&(Country[z]==Ctr)]))))
GLYuhum <- sort(unique(GLY[(Country[z]==Ctr)&(Group[z]=="Human")])) 
Nhumnovel[4] <- length(setdiff(GLYuhum,sort(unique(GLY[((Group[z]!="Human"))&(Country[z]==Ctr)]))))
PGMuhum <- sort(unique(PGM[(Country[z]==Ctr)&(Group[z]=="Human")]))
Nhumnovel[5] <- length(setdiff(PGMuhum,sort(unique(PGM[((Group[z]!="Human"))&(Country[z]==Ctr)]))))
TKTuhum <- sort(unique(TKT[(Country[z]==Ctr)&(Group[z]=="Human")])) 
Nhumnovel[6] <- length(setdiff(TKTuhum,sort(unique(TKT[((Group[z]!="Human"))&(Country[z]==Ctr)]))))
UNCuhum <- sort(unique(UNC[(Country[z]==Ctr)&(Group[z]=="Human")])) 
Nhumnovel[7] <- length(setdiff(UNCuhum,sort(unique(UNC[((Group[z]!="Human"))&(Country[z]==Ctr)]))))

# Nhum: number of all human isolates used as sample size in the binomial for number of mutations (novel types):
# for each locus:
for(i in 1:7){ Nhum[i] <- sum( (Country[z]==Ctr)&(Group[z]=="Human")  )   }


library("R2OpenBUGS")

beta <- rep(1,ns)  # Dir prior parameters for sources
alpha <- structure(.Data=rep(1,ns*ns),.Dim=c(ns,ns)) # Dir prior parameters for migration probabilities
data <- list("humansASP","humansGLN","humansGLT","humansGLY","humansPGM","humansTKT","humansUNC",
"sourcesASP","sourcesGLN","sourcesGLT","sourcesGLY","sourcesPGM","sourcesTKT","sourcesUNC",
"ns","nat","beta","alpha","Nhum","Nhumnovel",
"NASP","NGLN","NGLT","NGLY","NPGM","NTKT","NUNC",
"NASPnovel","NGLNnovel","NGLTnovel","NGLYnovel","NPGMnovel","NTKTnovel","NUNCnovel")

parameters <- c("qASP","phi","etaASP")
inits <- function(){list(g0=rep(1,ns),
pmuta=structure(.Data=rep(0.5,ns*7),.Dim=c(ns,7)),
h0=structure(.Data=rep(1,ns*ns),.Dim=c(ns,ns)))}

res2 <- bugs(data,inits,parameters,"SA_allelemodel2.txt",n.chains=1,n.burnin=100,n.iter=5100)
attach.bugs(res2)

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

parameters <- c("qASP","phi","etaASP")
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

#############
beta <- rep(1,ns)  # Dir prior parameters for sources
alpha <- structure(.Data=rep(1,ns*ns),.Dim=c(ns,ns)) # Dir prior parameters for migration probabilities
data <- list("humansASP","humansGLN","humansGLT","humansGLY","humansPGM","humansTKT","humansUNC",
"sourcesASP","sourcesGLN","sourcesGLT","sourcesGLY","sourcesPGM","sourcesTKT","sourcesUNC",
"ns","nat","beta","alpha","Nhum","Nhumnovel",
"NASP","NGLN","NGLT","NGLY","NPGM","NTKT","NUNC",
"NASPnovel","NGLNnovel","NGLTnovel","NGLYnovel","NPGMnovel","NTKTnovel","NUNCnovel")

parameters <- c("phi","etaASP")
inits <- function(){list(g0=rep(1,ns),
pmuta=structure(.Data=rep(0.5,ns*7),.Dim=c(ns,7)),
h0=structure(.Data=rep(1,ns*ns),.Dim=c(ns,ns)))}

res2B <- bugs(data,inits,parameters,"SA_allelemodel2B.txt",n.chains=1,n.burnin=100,n.iter=5100)
attach.bugs(res2B)







