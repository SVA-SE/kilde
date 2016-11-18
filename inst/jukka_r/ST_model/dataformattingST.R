
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

maxns <- ns+1  # one 'unknown source' added

# sample frequencies (counts) of each ST-type in humans and in sources:
sourcesST  <- matrix(0,maxns,length(STu)) 
humansST <- numeric()
for(i in 1:length(STu)){
humansST[i] <- sum((ST==STu[i])&(Group=="Human")&(Country==Ctr)&z)
for(j in 1:ns){
sourcesST[j,i] <- sum((ST==STu[i])&(Group==sourcenames[j])&(Country==Ctr)&z) }
}

# Number of all human isolates:
Nisolates <- sum((Country==Ctr)&(z)&(Group=="Human"))

# ST numbers for each human isolate:
HumanST <- ST[(Country==Ctr)&(z)&(Group=="Human")]

# Table of type indicators for isolates:
positionST <- 1:length(STu)
IST <- matrix(0,Nisolates,length(STu))
for(i in 1:Nisolates){
IST[i,positionST[STu==HumanST[i]]] <- 1 
}

# build index for ST types of the ith isolate:
ind <- numeric() 
prodIST<-numeric()
for(i in 1:Nisolates){
for(j in 1:length(STu)){prodIST[j] <- IST[i,j]*j}
ind[i] <- sum(prodIST[])
}

ns<-maxns  # one "unknown source" added, with data sample zero, or the following:


# Leave "unknown source" with zero observed counts,
# or set "data" for the unknown source as the average sample over all sources,
# or as the counts drawn from Unique Human Types:

if(UM==1){
for(i in 1:length(STu)){
 sourcesST[ns,i] <- round(mean(sourcesST[1:(ns-1),i]))  
}
}
if(UM==2){
for(i in 1:length(STu)){
sourcesST[ns,i] <- sum((ST==STu[i])&is.element(STu[i],STuHo)&(Group=="Human")&(Country==Ctr)&z) 
}
}


