library(kilde)
object <- clean_data(sample_data(), country = "Peru")
model <- system.file(package = "kilde", "june15/SA_allele_isolate.txt")
parameters <- c("qASP","phi","P","Z")

inits <- function(){list(g0=rep(1,ns))}
res <- run_model(model = model,
                 df = object,
                 inits = inits,
                 parameters = parameters)
)
## This must run clean
