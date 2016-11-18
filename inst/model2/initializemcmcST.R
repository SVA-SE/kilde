
# VARIABLE STRUCTURES FOR THE MCMC SAMPLER:
STpar0 <- matrix(0,ns,length(STu))
qST<- structure(.Data=c(rep(0,MCMC*ns*length(STu))),.Dim=c(MCMC,ns,length(STu)))
###qST <- matrix(NA,ns,length(STu))  # could be used if not saving the MCMC output
phii <- matrix(0,Nisolates,ns)
phii0 <- matrix(0,Nisolates,ns)
Z <- matrix(0,MCMC,Nisolates)
###Z <- numeric() # could be used if not saving the MCMC output
S <- matrix(0,Nisolates,ns) 
phi0 <- numeric()
phi <- matrix(0,MCMC,ns) 

###################
# INITIAL VALUES FOR MCMC:
for(i in 1:ns){
for(j in 1:length(STu)){qST[1,i,j] <- 1/length(STu)}
}
ga <- rgamma(ns,1,1)
phi[1,1:ns] <- ga/sum(ga) # source probabilities
Z[1,1:Nisolates]<-round(runif(Nisolates,0.5/ns,1+0.5/ns)*ns)
for(s in 1:Nisolates){
  for(k in 1:ns){
  S[s,k] <- Z[1,s]==k   # membership of sth isolate to group k  
  }
}

