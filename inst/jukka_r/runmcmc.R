##############################
##############################
## START MCMC GIBBS SAMPLER 
## *adapted * from 
## Pella & Masuda: Bayesian methods for analysis of stock
## mixtures from genetic characters. Fish. Bull. 99: 151-167 (2001).
h<-0  # control parameter (for avoiding very small values of q ~Dirich, maybe never needed)
## This code needs a series of objects to function. These need to be
## submitted to the code and are currently assumed to be in the global
## environment; They are:
## h, MCMC, ns, nat, FULL,
##
## qASP, ASPpar0, sourcesASP, IASP,
## qGLN, GLNpar0, sourcesGLN, IGLN,
## qGLT, GLTpar0, sourcesGLT, IGLT,
## qPGM, PGMpar0, sourcesPGM, IPGM,
## qGLY, GLYpar0, sourcesGLY, IGLY,
## qTKT, TKTpar0, sourcesTKT, ITKT,
## qUNC, UNCpar0, sourcesUNC, IUNC,
##
## Nisolates, phii0, ind, phii, S
## phi, phi0
##
## The follwoing of objects come from the data cleaning script:
##
## FULL, ind, and the sources*, and I*

for(mc in 2:MCMC){ 

    ## FULL CONDITIONAL DISTRIBUTIONS FOR 
    ## THE RELATIVE ALLELE TYPE FREQUENCIES IN EACH SOURCE POPULATION:
    ## (qASP, GLN, GLT, GLY, PGM, TKT, UNC: these are defined for each loci)
    for(i in 1:ns){   # ns = number of source populations
        for(j in 1:nat[1]){ 
            ASPpar0[i, j] <- rgamma(1,
                                    sourcesASP[i, j] + FULL * sum(S[, i] * IASP[, j] > 0) + 1/nat[1],
                                    1)
        }
        for(j in 1:nat[1]){
            qASP[mc, i, j] <-  (ASPpar0[i, j] + h) / (sum(ASPpar0[i, ]) + h * nat[1])
        }
        for(j in 1:nat[2]){ 
            GLNpar0[i, j] <- rgamma(1,
                                    sourcesGLN[i, j] + FULL * sum(S[, i] * IGLN[, j] > 0) + 1/nat[2],
                                    1)
        }
        for(j in 1:nat[2]){
            qGLN[mc, i, j] <-  (GLNpar0[i, j] + h) / (sum(GLNpar0[i, ]) + h * nat[2])
        }
        for(j in 1:nat[3]){ 
            GLTpar0[i, j] <- rgamma(1,
                                    sourcesGLT[i, j] + FULL * sum(S[, i] * IGLT[, j] > 0) + 1/nat[3],
                                    1)
        }
        for(j in 1:nat[3]){
            qGLT[mc, i, j] <-  (GLTpar0[i, j] + h) / (sum(GLTpar0[i, ]) + h * nat[3])
        }
        for(j in 1:nat[4]){ 
            GLYpar0[i, j] <- rgamma(1,
                                    sourcesGLY[i, j] + FULL * sum(S[, i] * IGLY[, j] > 0) + 1/nat[4],
                                    1)
        }
        for(j in 1:nat[4]){
            qGLY[mc, i, j] <-  (GLYpar0[i, j] + h) / (sum(GLYpar0[i, ]) + h * nat[4])
        }
        for(j in 1:nat[5]){ 
            PGMpar0[i, j] <- rgamma(1,
                                    sourcesPGM[i, j] + FULL * sum(S[, i] * IPGM[, j] > 0) + 1/nat[5],
                                    1)
        }
        for(j in 1:nat[5]){
            qPGM[mc, i, j] <-  (PGMpar0[i, j] + h) / (sum(PGMpar0[i, ]) + h * nat[5])
        }
        for(j in 1:nat[6]){ 
            TKTpar0[i, j] <- rgamma(1,
                                    sourcesTKT[i, j] + FULL * sum(S[, i] * ITKT[, j] > 0) + 1/nat[6],
                                    1)
        }
        for(j in 1:nat[6]){
            qTKT[mc, i, j] <-  (TKTpar0[i, j] + h) / (sum(TKTpar0[i, ]) + h * nat[6])
        }
        for(j in 1:nat[7]){ 
            UNCpar0[i, j] <- rgamma(1,
                                    sourcesUNC[i, j] + FULL * sum(S[, i] * IUNC[, j] > 0) + 1/nat[7],
                                    1)
        }
        for(j in 1:nat[7]){
            qUNC[mc, i, j] <-  (UNCpar0[i, j] + h) / (sum(UNCpar0[i, ]) + h * nat[7])
        }
    }  ## end of ns

    ## FULL CONDITIONAL DISTRIBUTION FOR EACH Z: 
    ## (Z = group indicators for each isolate) 
    for(i in 1:Nisolates){
        for(k in 1:ns){ 
            phii0[i, k] <- exp( log(phi[mc-1, k]) +
                                log(qASP[mc, k, ind[i, 1]]) +
                                log(qGLN[mc, k, ind[i, 2]]) +
                                log(qGLT[mc, k, ind[i, 3]]) +
                                log(qGLY[mc, k, ind[i, 4]]) +
                                log(qPGM[mc, k, ind[i, 5]]) +
                                log(qTKT[mc, k, ind[i, 6]]) +
                                log(qUNC[mc, k, ind[i, 7]])
                               )
        }
        for(k in 1:ns){
            phii[i, k] <- exp(log(phii0[i, k])-log(sum(phii0[i, ])) ) ## normalize
        }
        nn <- 1:ns
        Z[mc, i] <- min(nn[runif(1)<cumsum(phii[i, ])])     ## "rcat(1, phii[i, ])"
        
        for(k in 1:ns){
            S[i, k] <- Z[mc, i]==k   ## membership of ith isolate to group k  
        }
    } ## end Nisolates

    ## FULL CONDITIONAL DISTRIBUTION FOR phi:
    ## (phi = proportions of source populations)
    for(k in 1:ns){
        phi0[k] <- rgamma(1,
                          sum(S[, k]) + 1/ns,
                          1)
    }
    for(k in 1:ns){
        phi[mc, k] <- phi0[k]/sum(phi0[])
    }  
}  ## end of MCMC
################################################
