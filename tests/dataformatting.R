library(kilde)
df <- sample_data()
df <- df[df$country == "Canada",]
load(system.file("extdata/dataformatting_ob.Rda", package = "kilde"))
ob1 <- dataformatting(DATA = df,
                     UM = 2)
stopifnot(identical(ob, ob1))
rm(list = ls())
