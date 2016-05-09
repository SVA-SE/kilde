# A package to perform Bayesian source attribution of Campylobacter in R

An R package that connects to
[JAGS](http://mcmc-jags.sourceforge.net/) via
[rjags](https://cran.r-project.org/web/packages/rjags/index.html) to
perform source attribution analysis on MLST data from *Campylobacter
jejuni*


## Installation

```
library(devtools)
install_github("trosendal/kilde")
```

## Example

```
library(kilde)
df <- sample_data()
```
