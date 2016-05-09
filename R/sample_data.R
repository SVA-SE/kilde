##' Read in the sample data and clean it
##' 
##' @title sample_data 
##' @return A data.frame
##' @author Thomas Rosendal
##' @export
sample_data <- function(){
    df <- read.table(system.file("extdata/sampledata.txt", package = "kilde"),
                     header = TRUE)
    return(df)
}
