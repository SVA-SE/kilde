##############################
##############################
# START MCMC 
# GIBBS SAMPLER *adapted* from 
# Pella & Masuda: Bayesian methods for analysis of stock
# mixtures from genetic characters. Fish. Bull. 99: 151-167 (2001).
h<-0  # control parameter (for avoiding very small values of q ~Dirich, maybe not needed)
for(mc in 2:MCMC){ 

# FULL CONDITIONAL DISTRIBUTIONS FOR 
# THE RELATIVE ST TYPE FREQUENCIES IN EACH SOURCE POPULATION:
for(i in 1:ns){   # ns = number of source populations
for(j in 1:length(STu)){ 
  STpar0[i,j] <- rgamma(1,sourcesST[i,j]+FULL*sum(S[,i]*IST[,j]>0)+1/length(STu),1)
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
Z[mc,i] <- min(nn[runif(1)<cumsum(phii[i,])])     # "rcat(1,phii[i,])"
 
for(k in 1:ns){
    S[i,k] <- Z[mc,i]==k   # membership of ith isolate to group k  
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
