#########
#### ST-model

setwd("C://Users/Jukkis/Documents/Evira/Kampy/SourceAttribution/Mallit/Lopulliset/ToThomas")
DATA <- read.table("NMDD2015data_vers21_modif.txt",header=TRUE)

attach(DATA)
# use ALL allele types
z <- !is.na(ST)

Ctr = "SE"
#sourcenames=c("Cattle","Broiler","Turkey","Water") # FI  
sourcenames=c("Bird","Cattle","Sheep","Broiler","Pig","Dog","Water") # SE 
#sourcenames=c("Cattle","Broiler","Pig","Turkey","Duck") # DK
#sourcenames=c("Broiler","Sheep","Pig","Cattle") # NO
ns <- length(sourcenames)

# find how many different ST-types there are (in all source groups), and list them:
STinhumans <- sort(unique(ST[(Country==Ctr)&(z)&(Group=="Human")]))
STinsources <- sort(unique(ST[(Country==Ctr)&(z)&(Group!="Human")]))
STu <- STinsources


# find how many different ST-types there are (only in humans), and list them:
STuHo <- setdiff(sort(unique(ST[(Country==Ctr)&(z)&(Group=="Human")])),STinsources)


# number of all different STs (in all sources)
nST <- length(STu) 

# number of all different STs (only found in humans)
nSTHo <- length(STuHo) 

# number (counts) of ST-types found in sources:
sourcesST <- matrix(NA,ns,length(STu)) 
# number (counts) of such ST-types in humans:
humansST <- numeric()

# sample frequencies (counts) of each ST-type (found in sources)
for(i in 1:length(STu)){
humansST[i] <- sum((ST==STu[i])&(Group=="Human")&(Country==Ctr)&z)
for(j in 1:ns){
sourcesST[j,i] <- sum((ST==STu[i])&(Group==sourcenames[j])&(Country==Ctr)&z) }
}


# number (counts) of ST-types found ONLY in humans:
humansSTo <- numeric()

# sample frequencies (counts) of each ST-type (the ST-types found only in humans):
for(i in 1:length(STuHo)){
humansSTo[i] <- sum((ST==STuHo[i])&(Group=="Human")&(Country==Ctr)&z)
}

# number of all human isolates:
Nhumans <- sum((Country==Ctr)&(z)&(Group=="Human"))

# number of human isolates with an ST-type not found in sources:
# (out of a total number of Nhumans). 
# To be used in BUGS for estimating proportion of "unknown sources"
NSTHo <- sum(humansSTo)

# "NST" is to be used as sample size in the binomial model for 
# "mutations" ("novel" ST-type or "distinct" ST-type: not found in other sources) 
# per source (for each locus):
NST <- numeric()

for(i in 1:ns){
NST[i] <- sum( (Country==Ctr)&(Group==sourcenames[i])&z  ) 
}

# BELOW: loop over sources, list ST types in each source

maxlength <- numeric(); STle <-numeric()

for(i in 1:ns){ 
 STle[i]<-length(sort(unique(ST[(Group==sourcenames[i])&(Country==Ctr)&z])))
 }
maxlength <- max(STle)
STusource <- matrix(NA,ns,maxlength)

# the list of found ST-types for each source:
# STusource[] contains listing of ST-types found in each source
for(i in 1:ns){  
STusource[i,1:STle[i]] <- sort(unique(ST[(Group==sourcenames[i])&(Country==Ctr)&z])) 
}

STnovel <- matrix(NA,ns,maxlength)
STnn <- numeric()
  
# Find out list of 'novel types' in a source (that were not found in other sources)
for(i in 1:ns){ 
STnn[i]<-length(setdiff(STusource[i,1:STle[i]],sort(unique(ST[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z]))))
if(STnn[i]>0){
STnovel[i,1:STnn[i]] <- setdiff(STusource[i,1:STle[i]],sort(unique(ST[((Group!=sourcenames[i])&(Group!="Human"))&(Country==Ctr)&z])))
}
}

# count the frequencies (numbers) of the novel types:
sourcesSTnovel=matrix(0,ns,nST)

for(j in 1:ns){
 for(i in 1:nST){
 sourcesSTnovel[j,i] <- sourcesST[j,i]*is.element(STu[i],STnovel[j,1:STnn[j]]) }
 }   

# number of all 'novel types' per source for BUGS-model:  
# (ST-types that were not in other sources)
NSTnovel <- numeric()
for(i in 1:ns){ NSTnovel[i] <- sum(sourcesSTnovel[i,]) }

### find out what ST-types there are in a source that were found at least 
# in some other sources:
STnco <- numeric()   
for(i in 1:ns){
STnco[i] <- length(setdiff(STusource[i,],STnovel[i,]))
}
STcommon <- matrix(NA,ns,max(STnco))

for(i in 1:ns){
STcommon[i,1:STnco[i]] <- setdiff(STusource[i,],STnovel[i,])
}
migrateST <- matrix(NA,ns,nST)

# count frequencies of "migrated" types per source:
for(j in 1:ns){
for(i in 1:nST){
migrateST[j,i] <- sourcesST[j,i]*is.element(STu[i],STcommon[j,])
}
}

library("R2OpenBUGS")

# the migration-mutation-unknownsource -model:

beta <- rep(1,ns) # Dir prior parameters for sources
alpha <- structure(.Data=rep(1,ns*ns),.Dim=c(ns,ns)) # Dir prior parameters for migration probabilities
data <- list("humansST",
"migrateST",
"sourcesSTnovel",
"sourcesST","ns","nST","beta","alpha",
"NST","NSTnovel",
"NSTHo","Nhumans")

parameters <- c("phi","etaST","khi","pmuta","pother")
inits <- function(){list(g0=rep(1,ns),
pmuta=rep(0.5,ns),
h0=structure(.Data=rep(1,ns*ns),.Dim=c(ns,ns)))}

resnewNO <- bugs(data,inits,parameters,"SA_STmodel2New.txt",n.chains=1,n.burnin=100,n.iter=5100)
attach.bugs(resnewNO)


mm <- numeric();for(i in 1:ns){mm[i]<-mean(phi[,i]*(1-pother))}
mmo <- order(mm)
ch <- mmo[ns]

par(mfrow=c(ns,1))
for(i in 1:ns){
plot(density(phi[,i]*(1-pother)),main="marginal posterior",xlab=sourcenames[i],xlim=c(0,1));points(mean(phi[,i]*(1-pother)),0,col="red",pch=16) 
if(i==ch){points(mean(phi[,i]*(1-pother)),0,col="blue",pch=16,cex=2);text(mean(phi[,i]*(1-pother)),1,"MaxP",cex=1.3)}
}







