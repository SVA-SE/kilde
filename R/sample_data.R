##' Read in the sample data and clean it
##' 
##' @title sample_data 
##' @return A data.frame
##' @author Thomas Rosendal
##' @export
sample_data <- function(){
    df <- read.table(system.file("extdata/sampledata.txt", package = "kilde"),
                     header = TRUE)
    df$group <- factor(df$group,
                       levels = 0:9,
                       labels = c("human",
                                  "alligator",
                                  "baboon",
                                  "whale",
                                  "parrot",
                                  "panda",
                                  "frog",
                                  "rabbit",
                                  "pig",
                                  "lion"))
    df$country <- factor(df$country,
                         levels = 1:4,
                         labels = c("Peru",
                                    "Canada",
                                    "Estonia",
                                    "Japan"))
    return(df)
}
