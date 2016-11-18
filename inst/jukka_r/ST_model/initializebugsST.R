
# BUGS-model with multinomial(q[Z],1)-likelihoods 
# for each human isolate with unknown group indicator (Z) 

nt <- length(STu)
sources <- sourcesST
Nisolates <- sum(humansST)
beta <- rep(1/ns,ns) # Dir prior parameters for sources
data <- list("IST","sourcesST","ns","nt","beta","Nisolates","FULL")

parameters <- c("phi","Z","qST")
inits <- function(){list(g0=rgamma(ns,1,1),
q0=structure(.Data=rep(1,ns*nt),.Dim=c(ns,nt)),
Z=round(runif(Nisolates,0.5/ns,1+0.5/ns)*ns))}
