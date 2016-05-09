##' Clean data to be passed on to the model function
##' 
##' @title clean_data 
##' @param df 
##' @return list
##' @author Thomas Rosendal
##' @export
clean_data <- function(df = sample_data(), country){
    df <- df[df$country == country,]
    ASP <- table(df$group, df$ASP)
    GLN <- table(df$group, df$GLN)
    GLT <- table(df$group, df$GLT)
    GLY <- table(df$group, df$GLY)
    PGM <- table(df$group, df$PGM)
    TKT <- table(df$group, df$TKT)
    UNC <- table(df$group, df$UNC)
    ns <- length(unique(df$group[df$group != 0]))
    beta <- rep(1,ns)
    nat <- c(length(unique(df$ASP)),
             length(unique(df$GLN)),
             length(unique(df$GLT)),
             length(unique(df$GLY)),
             length(unique(df$PGM)),
             length(unique(df$TKT)),
             length(unique(df$UNC)))
    alpha <- structure(.Data=rep(1,ns*ns),.Dim=c(ns,ns))
    result <- list(ASP = ASP, GLN = GLN, GLT = GLT, GLY = GLY,
                   PGM = PGM, TKT = TKT, UNC = UNC, ns =  ns,
                   beta = beta, nat = nat, alpha = alpha)
## The follwing are arbitrary changes to fulfill the structure of the
## BUGS code
    humansASP <- as.vector(result$ASP[1,])
    humansGLN <- as.vector(result$GLN[1,])
    humansGLT <- as.vector(result$GLT[1,])
    humansGLY <- as.vector(result$GLY[1,])
    humansPGM <- as.vector(result$PGM[1,])
    humansTKT <- as.vector(result$TKT[1,])
    humansUNC <- as.vector(result$UNC[1,])
    sourcesASP <- as.matrix(result$ASP[2:(ns+1),])
    sourcesGLN <- as.matrix(result$GLN[2:(ns+1),])
    sourcesGLT <- as.matrix(result$GLT[2:(ns+1),])
    sourcesGLY <- as.matrix(result$GLY[2:(ns+1),])
    sourcesPGM <- as.matrix(result$PGM[2:(ns+1),])
    sourcesTKT <- as.matrix(result$TKT[2:(ns+1),])
    sourcesUNC <- as.matrix(result$UNC[2:(ns+1),])
    result <- list(humansASP = humansASP,
                   humansGLN = humansGLN,
                   humansGLT = humansGLT,
                   humansGLY = humansGLY,
                   humansPGM = humansPGM,
                   humansTKT = humansTKT,
                   humansUNC = humansUNC,
                   sourcesASP = sourcesASP,
                   sourcesGLN = sourcesGLN,
                   sourcesGLT = sourcesGLT,
                   sourcesGLY = sourcesGLY,
                   sourcesPGM = sourcesPGM,
                   sourcesTKT = sourcesTKT,
                   sourcesUNC = sourcesUNC,
                   ns = ns,
                   nat = nat,
                   beta = beta,
                   alpha = alpha)
    
    return(result)
}
