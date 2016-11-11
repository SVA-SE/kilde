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
## The following of objects come from the data cleaning script:
##
## FULL, ind, and the sources*, and I*

runmcmc <- function(result,
                    ob,
                    MCMC,
                    h = 0,
                    FULL = 0){
    
    for(mc in 2:MCMC){ 

        ## FULL CONDITIONAL DISTRIBUTIONS FOR 
        ## THE RELATIVE ALLELE TYPE FREQUENCIES IN EACH SOURCE POPULATION:
        ## (qASP, GLN, GLT, GLY, PGM, TKT, UNC: these are defined for each loci)
        for(i in 1:ob$inits$ns){   # ns = number of source populations
            for(j in 1:ob$inits$nat[1]){ 
                result$ASPpar0[i, j] <- rgamma(1,
                                        ob$data$sourcesASP[i, j] + FULL * sum(result$S[, i] * ob$data$IASP[, j] > 0) + 1/ob$inits$nat[1],
                                        1)
            }
            for(j in 1:ob$inits$nat[1]){
                result$qASP[mc, i, j] <-  (result$ASPpar0[i, j] + h) / (sum(result$ASPpar0[i, ]) + h * ob$inits$nat[1])
            }
            for(j in 1:ob$inits$nat[2]){ 
                result$GLNpar0[i, j] <- rgamma(1,
                                        ob$data$sourcesGLN[i, j] + FULL * sum(result$S[, i] * ob$data$IGLN[, j] > 0) + 1/ob$inits$nat[2],
                                        1)
            }
            for(j in 1:ob$inits$nat[2]){
                result$qGLN[mc, i, j] <-  (result$GLNpar0[i, j] + h) / (sum(result$GLNpar0[i, ]) + h * ob$inits$nat[2])
            }
            for(j in 1:ob$inits$nat[3]){ 
                result$GLTpar0[i, j] <- rgamma(1,
                                        ob$data$sourcesGLT[i, j] + FULL * sum(result$S[, i] * ob$data$IGLT[, j] > 0) + 1/ob$inits$nat[3],
                                        1)
            }
            for(j in 1:ob$inits$nat[3]){
                result$qGLT[mc, i, j] <-  (result$GLTpar0[i, j] + h) / (sum(result$GLTpar0[i, ]) + h * ob$inits$nat[3])
            }
            for(j in 1:ob$inits$nat[4]){ 
                result$GLYpar0[i, j] <- rgamma(1,
                                        ob$data$sourcesGLY[i, j] + FULL * sum(result$S[, i] * ob$data$IGLY[, j] > 0) + 1/ob$inits$nat[4],
                                        1)
            }
            for(j in 1:ob$inits$nat[4]){
                result$qGLY[mc, i, j] <-  (result$GLYpar0[i, j] + h) / (sum(result$GLYpar0[i, ]) + h * ob$inits$nat[4])
            }
            for(j in 1:ob$inits$nat[5]){ 
                result$PGMpar0[i, j] <- rgamma(1,
                                        ob$data$sourcesPGM[i, j] + FULL * sum(result$S[, i] * ob$data$IPGM[, j] > 0) + 1/ob$inits$nat[5],
                                        1)
            }
            for(j in 1:ob$inits$nat[5]){
                result$qPGM[mc, i, j] <-  (result$PGMpar0[i, j] + h) / (sum(result$PGMpar0[i, ]) + h * ob$inits$nat[5])
            }
            for(j in 1:ob$inits$nat[6]){ 
                result$TKTpar0[i, j] <- rgamma(1,
                                        ob$data$sourcesTKT[i, j] + FULL * sum(result$S[, i] * ob$data$ITKT[, j] > 0) + 1/ob$inits$nat[6],
                                        1)
            }
            for(j in 1:ob$inits$nat[6]){
                result$qTKT[mc, i, j] <-  (result$TKTpar0[i, j] + h) / (sum(result$TKTpar0[i, ]) + h * ob$inits$nat[6])
            }
            for(j in 1:ob$inits$nat[7]){ 
                result$UNCpar0[i, j] <- rgamma(1,
                                        ob$data$sourcesUNC[i, j] + FULL * sum(result$S[, i] * ob$data$IUNC[, j] > 0) + 1/ob$inits$nat[7],
                                        1)
            }
            for(j in 1:ob$inits$nat[7]){
                result$qUNC[mc, i, j] <-  (result$UNCpar0[i, j] + h) / (sum(result$UNCpar0[i, ]) + h * ob$inits$nat[7])
            }
        }  ## end of ns

        ## FULL CONDITIONAL DISTRIBUTION FOR EACH Z: 
        ## (Z = group indicators for each isolate) 
        for(i in 1:result$Nisolates){
            for(k in 1:ob$inits$ns){ 
                result$phii0[i, k] <- exp( log(result$phi[mc-1, k]) +
                                    log(result$qASP[mc, k, ob$data$ind[i, 1]]) +
                                    log(result$qGLN[mc, k, ob$data$ind[i, 2]]) +
                                    log(result$qGLT[mc, k, ob$data$ind[i, 3]]) +
                                    log(result$qGLY[mc, k, ob$data$ind[i, 4]]) +
                                    log(result$qPGM[mc, k, ob$data$ind[i, 5]]) +
                                    log(result$qTKT[mc, k, ob$data$ind[i, 6]]) +
                                    log(result$qUNC[mc, k, ob$data$ind[i, 7]])
                                   )
            }
            for(k in 1:ob$inits$ns){
                result$phii[i, k] <- exp(log(result$phii0[i, k])-log(sum(result$phii0[i, ])) ) ## normalize
            }
            nn <- 1:ob$inits$ns
            result$Z[mc, i] <- min(nn[runif(1) < cumsum(result$phii[i, ])])     ## "rcat(1, result$phii[i, ])"
            
            for(k in 1:ob$inits$ns){
                result$S[i, k] <- result$Z[mc, i]==k   ## membership of ith isolate to group k  
            }
        } ## end Nisolates

        ## FULL CONDITIONAL DISTRIBUTION FOR phi:
        ## (phi = proportions of source populations)
        for(k in 1:ob$inits$ns){
            result$phi0[k] <- rgamma(1,
                              sum(result$S[, k]) + 1/ob$inits$ns,
                              1)
        }
        for(k in 1:ob$inits$ns){
            result$phi[mc, k] <- result$phi0[k]/sum(result$phi0[])
        }  
    }  ## end of MCMC
    result <- list(var_a = result, var_b = ob)
    return(result)
}
