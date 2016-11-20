##' Data formatting for the ST based model
##'
##' @title data_formatting
##' @param DATA A data.frame with at least A "Group" and "ST" column
##' @param UM An integer. Choose 0 for setting 0 as 'data sample' for
##'     unknown source.  Choose 1 for taking the mean type frequencies
##'     as 'data sample' for unknown source.  Choose 2 for taking the
##'     type frequencies drawn from the STs that were unique to
##'     humans, as 'data sample' for unknown source.
##' @return A list that may be submitted to the function to initialize
##'     the models in R or bugs
##' @author Jukka Ranta
##' @export
dataformatting_ST <- function(DATA, UM = 2) {
    z <- !is.na(DATA$ST)
    sourcenames <- setdiff(unique(DATA$group), "human")
    ns <- length(sourcenames)

    ## find how many different STs there are (in all isolates), and
    ## list them:
    STu <- sort(unique(DATA$ST[z]))
    STuH <- sort(unique(DATA$ST[z & (DATA$group == "human")]))
    STuS <- sort(unique(DATA$ST[z & (DATA$group != "human")]))
    STuHo <- setdiff(STuH, STuS) 

    HumanST <- DATA$ST[z & (DATA$group == "human")]
    Humnovel <- 0
    for(i in 1:length(HumanST)){
        Humnovel <- Humnovel + is.element(HumanST[i] , STuHo)*1
    }
    ## Sample probability to see ST in Human isolates that was not
    ## seen in sources:
    PUNST <- Humnovel / length(HumanST)

    maxns <- ns + 1  ## one 'unknown source' added

    ## sample frequencies (counts) of each ST-type in humans and in
    ## sources:
    sourcesST  <- matrix(0, maxns, length(STu)) 
    humansST <- numeric()
    for(i in 1:length(STu)){
        humansST[i] <- sum((DATA$ST == STu[i]) & (DATA$group == "human") & z)
        for(j in 1:ns){
            sourcesST[j,i] <- sum((DATA$ST == STu[i]) & (DATA$group == sourcenames[j]) & z)
        }
    }

    ## Number of all human isolates:
    Nisolates <- sum(z & (DATA$group == "human"))

    ## ST numbers for each human isolate:
    HumanST <- DATA$ST[z & (DATA$group == "human")]

    ## Table of type indicators for isolates:
    positionST <- 1:length(STu)
    IST <- matrix(0, Nisolates, length(STu))
    for(i in 1:Nisolates){
        IST[i, positionST[STu == HumanST[i]]] <- 1 
    }

    ## build index for ST types of the ith isolate:
    ind <- numeric() 
    prodIST<-numeric()
    for(i in 1:Nisolates){
        for(j in 1:length(STu)){
            prodIST[j] <- IST[i, j] * j
        }
        ind[i] <- sum(prodIST[])
    }

    ns <- maxns  ## one "unknown source" added, with data sample zero,
    ## or the following:

    ## Leave "unknown source" with zero observed counts, or set "data"
    ## for the unknown source as the average sample over all sources,
    ## or as the counts drawn from Unique Human Types:

    if(UM == 1){
        for(i in 1:length(STu)){
            sourcesST[ns, i] <- round(mean(sourcesST[1:(ns - 1), i]))  
        }
    }
    if(UM == 2){
        for(i in 1:length(STu)){
            sourcesST[ns, i] <- sum((DATA$ST == STu[i]) &
                                    is.element(STu[i], STuHo) &
                                    (DATA$group == "human") &
                                    z)
        }
    }
    result <- list(DATA = DATA,
                   humansST = humansST,
                   HumanST = HumanST,
                   Humnovel = Humnovel,
                   ind = ind,
                   IST = IST,
                   maxns = maxns,
                   Nisolates = Nisolates,
                   ns = ns,
                   positionST = positionST,
                   prodIST = prodIST,
                   PUNST = PUNST,
                   sourcenames = sourcenames,
                   sourcesST = sourcesST,
                   STu = STu,
                   STuH = STuH,
                   STuHo = STuHo,
                   STuS = STuS,
                   UM = UM)
    return(result)
}
