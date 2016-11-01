library(kilde)

## Test 0 - A clean run of a model

df <- clean_data(sample_data(), "Canada")
inits <- function() {
    list(g0 = rep(1, df$bugs_data$ns), 
         pmuta = structure(.Data = rep(0.5, df$bugs_data$ns * 7), 
                           .Dim = c(df$bugs_data$ns, 7)), 
         h0=structure(.Data=rep(1, df$bugs_data$ns * df$bugs_data$ns), 
                      .Dim=c(df$bugs_data$ns, df$bugs_data$ns)))
}

parameters <- c("qASP", "phi", "etaASP")

run_model(model = "SA_allelemodel1.jag",
          df,
          inits,
          parameters,
          n.chains = 2,
          n.burnin = 99,
          n.iter = 999)

## Test 1 - The data object has the wrong class
df <- clean_data(sample_data(), "Canada")
class(df) <- "foo"
inits <- function() {
    list(g0 = rep(1, df$bugs_data$ns), 
         pmuta = structure(.Data = rep(0.5, df$bugs_data$ns * 7), 
                           .Dim = c(df$bugs_data$ns, 7)), 
         h0=structure(.Data=rep(1, df$bugs_data$ns * df$bugs_data$ns), 
                      .Dim=c(df$bugs_data$ns, df$bugs_data$ns)))
}

parameters <- c("qASP", "phi", "etaASP")

res <- tools::assertError(
    run_model(model = "SA_allelemodel1.jag",
              df,
              inits,
              parameters)
)
stopifnot(length(grep("The class of the input object 'df' must be c\\('list','kilde_data'\\) which indicates that it comes from the clean_data\\(\\) function in the kilde package",
                     res[[1]]$message)) > 0)
rm(list = ls())

## Test 2 - The model file doesn't exist

df <- clean_data(sample_data(), "Canada")
inits <- function() {
    list(g0 = rep(1, df$bugs_data$ns), 
         pmuta = structure(.Data = rep(0.5, df$bugs_data$ns * 7), 
                           .Dim = c(df$bugs_data$ns, 7)), 
         h0=structure(.Data=rep(1, df$bugs_data$ns * df$bugs_data$ns), 
                      .Dim=c(df$bugs_data$ns, df$bugs_data$ns)))
}

parameters <- c("qASP", "phi", "etaASP")

res <- tools::assertError(
    run_model(model = "SA_allelemodelfoo.jag",
              df,
              inits,
              parameters)
)
stopifnot(length(grep("The model argument must be either the name of one of the models in the package or the complete path to a modelfile elsewhere",
                     res[[1]]$message)) > 0)
rm(list = ls())

## Test 3 - The inits is not a function

df <- clean_data(sample_data(), "Canada")
inits <- function() {
    list(g0 = rep(1, df$bugs_data$ns), 
         pmuta = structure(.Data = rep(0.5, df$bugs_data$ns * 7), 
                           .Dim = c(df$bugs_data$ns, 7)), 
         h0=structure(.Data=rep(1, df$bugs_data$ns * df$bugs_data$ns), 
                      .Dim=c(df$bugs_data$ns, df$bugs_data$ns)))
}
class(inits) <- "foo"

parameters <- c("qASP", "phi", "etaASP")

res <- tools::assertError(
    run_model(model = "SA_allelemodel.jag",
              df,
              inits,
              parameters)
)
stopifnot(length(grep("inits must be a function that returns the inits object",
                     res[[1]]$message)) > 0)
rm(list = ls())

## Test 4 - The chains burnin or iter are not numbers

df <- clean_data(sample_data(), "Canada")
inits <- function() {
    list(g0 = rep(1, df$bugs_data$ns), 
         pmuta = structure(.Data = rep(0.5, df$bugs_data$ns * 7), 
                           .Dim = c(df$bugs_data$ns, 7)), 
         h0=structure(.Data=rep(1, df$bugs_data$ns * df$bugs_data$ns), 
                      .Dim=c(df$bugs_data$ns, df$bugs_data$ns)))
}

parameters <- c("qASP", "phi", "etaASP")

res <- tools::assertError(
    run_model(model = "SA_allelemodel.jag",
              df,
              inits,
              parameters,
              n.chains = "s")
)
stopifnot(length(grep("n.chains must be numeric or integer",
                     res[[1]]$message)) > 0)

res <- tools::assertError(
    run_model(model = "SA_allelemodel.jag",
              df,
              inits,
              parameters,
              n.burnin = "s")
)
stopifnot(length(grep("n.burnin must be numeric or integer",
                     res[[1]]$message)) > 0)

res <- tools::assertError(
    run_model(model = "SA_allelemodel.jag",
              df,
              inits,
              parameters,
              n.iter = "s")
)
stopifnot(length(grep("n.iter must be numeric or integer",
                     res[[1]]$message)) > 0)

rm(list = ls())
