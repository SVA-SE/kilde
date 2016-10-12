library(kilde)

## Test 1 Read and clean the data in the expected way
df <- sample_data()
country <- "Canada"
clean_data(df, country)
rm(list = ls())

## Test 2 - REad and clean the data but provide the wrong columns
df <- sample_data()
names(df)[1] <- "foo"
country <- "Canada"
res <- tools::assertError(
    clean_data(df, country)
)
stopifnot(length(grep("The input data must have the column names:\'id', 'country', 'group', 'ST', 'ASP',
                               'GLN', 'GLT', 'GLY', 'PGM', 'TKT',
                               'UNC'",
                     res[[1]]$message)) > 0)
rm(list = ls())
## Test 3 - Read and clean the data but provide the wrong class of the group variable
df <- sample_data()
country <- "Canada"
class(df$group) <- "foo"
res <- tools::assertError(
    clean_data(df, country)
)
stopifnot(length(grep("The group variable should be a factor",
                     res[[1]]$message)) > 0)
rm(list = ls())
## Test 4 - Read and clean the data but provide the wrong class of the country variable
df <- sample_data()
country <- "Canada"
class(df$country) <- "foo"
res <- tools::assertError(
    clean_data(df, country)
)
stopifnot(length(grep("The country variable should be a factor or a character",
                     res[[1]]$message)) > 0)
rm(list = ls())
## Test 5 - select a country that has no data
df <- sample_data()
country <- "China"
res <- tools::assertError(
    clean_data(df, country)
)
stopifnot(length(grep("There are no observations in the data with the selected country name",
                     res[[1]]$message)) > 0)
rm(list = ls())
