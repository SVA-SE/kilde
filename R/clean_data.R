##' Clean data to be passed on to the model function
##' 
##' @title clean_data 
##' @param df 
##' @return list
##' @author Thomas Rosendal
##' @export
clean_data <- function(df, country){
    df <- df[df$country == country,]
    ASP <- table(df$group, df$ASP)
    GLN <- table(df$group, df$GLN)
    GLT <- table(df$group, df$GLT)
    GLY <- table(df$group, df$GLY)
    PGM <- table(df$group, df$PGM)
    TKT <- table(df$group, df$TKT)
    UNC <- table(df$group, df$UNC)
    ns <- length(unique(df$group))
    beta <- rep(1,ns)
    nat <- c(length(unique(df$ASP)),
             length(unique(df$GLN)),
             length(unique(df$GLT)),
             length(unique(df$GLY)),
             length(unique(df$PGM)),
             length(unique(df$TKT)),
             length(unique(df$UNC)))
    result <- list(ASP,
                   GLN,
                   GLT,
                   GLY,
                   PGM,
                   TKT,
                   UNC,
                   ns,
                   beta,
                   nat)
    return(result)
}
