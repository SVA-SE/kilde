## Copyright (c) 2016, Jukka Ranta, Thomas Rosendal
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
## 1. Redistributions of source code must retain the above copyright
## notice, this list of conditions and the following disclaimer.
##
## 2. Redistributions in binary form must reproduce the above copyright
## notice, this list of conditions and the following disclaimer in the
## documentation and/or other materials provided with the distribution.
##
## 3. Neither the name of the copyright holder nor the names of its
## contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##' Clean data to be passed on to the model function
##'
##' This function takes a dataset that has the Structure :
##'
##' id: A number or character that is the identifier of the individual
##' sample. This is actually not used by the function, but required as
##' it is passed onto the modelling tool.
##'
##' country: A character or factor of the county of origin that the
##' analysis should be compleated for. If the dataset only has one
##' country just supply a variable with one repeated value.
##'
##' group: Must be a factor. The first level of the factor must be the
##' reference group. This will frequently be humans.
##'
##' ST: An integer value repersenting the type identifier
##'
##' ASP, GLN, GLT, GLY, PGM, TKT, UNC: Integer values representing the
##' type identifier for each of the 7 MLST genes. 
##' 
##' @title clean_data 
##' @param df The dataframe containing the samples
##' @param country An country to subset the data by 
##' @return A list of which the first item goes to bugs
##' @author Thomas Rosendal
##' @export
clean_data <- function(df = sample_data(), country){

    if(!identical(names(df), c("id", "country", "group", "ST", "ASP",
                               "GLN", "GLT", "GLY", "PGM", "TKT",
                               "UNC"))) {
        stop("The input data must have the column names:\'id', 'country', 'group', 'ST', 'ASP',
                               'GLN', 'GLT', 'GLY', 'PGM', 'TKT',
                               'UNC'")
    }
    if(!identical(class(df$group), "factor")) {
        stop("The group variable should be a factor")
    }
    if(!(class(df$country) %in% c("factor", "character"))) {
        stop("The country variable should be a factor or a character")
    }

    if(!(identical(class(df$ST), "integer"))) {
        stop("ST must be an integer")
    }
    if(!(identical(class(df$ASP), "integer"))) {
        stop("ASP must be an integer")
    }
    if(!(identical(class(df$GLN), "integer"))) {
        stop("GLN must be an integer")
    }
    if(!(identical(class(df$GLT), "integer"))) {
        stop("GLT must be an integer")
    }
    if(!(identical(class(df$GLY), "integer"))) {
            stop("GLY must be an integer")
    }
    if(!(identical(class(df$PGM), "integer"))) {
        stop("PGM must be an integer")
    }
    if(!(identical(class(df$TKT), "integer"))) {
        stop("TKT must be an integer")
    }
    if(!(identical(class(df$UNC), "integer"))) {
        stop("UNC must be an integer")
    }

    df <- df[df$country == country,]

    if((nrow(df) == 0)) {
        stop("There are no observations in the data with the selected country name")
    }
    
    ASP <- table(df$group, df$ASP)
    GLN <- table(df$group, df$GLN)
    GLT <- table(df$group, df$GLT)
    GLY <- table(df$group, df$GLY)
    PGM <- table(df$group, df$PGM)
    TKT <- table(df$group, df$TKT)
    UNC <- table(df$group, df$UNC)
    ## Which groups have observations?
    group_names <- levels(df$group)[table(df$group) != 0]
    ## Which groups are not the first (The first should be humans or
    ## are at least considered the reference group)
    group_names <- group_names[2:length(group_names)]
    ns <- length(group_names)
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
                   beta = beta, nat = nat, alpha = alpha,
                   group_names = group_names)
    ## The following are arbitrary changes to fulfil the structure of
    ## the BUGS code
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
    result <- list(bugs_data = result,
                   group_names = group_names)
    class(result) <- c(class(result), "kilde_data")
    return(result)
}
