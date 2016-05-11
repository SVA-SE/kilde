|Linux Build|Windows Build|
|----|----|
|[![Build Status](https://travis-ci.org/trosendal/kilde.svg)](https://travis-ci.org/trosendal/kilde)|[![Build status](https://ci.appveyor.com/api/projects/status/dlxa7wxqs2lrc2eh?svg=true)](https://ci.appveyor.com/project/trosendal/kilde)|
# A package to perform Bayesian source attribution of Campylobacter in R

An R package that connects to
[OpenBUGS](http://www.openbugs.net) via
[R2OpenBUGS](https://cran.r-project.org/web/packages/R2OpenBUGS/vignettes/R2OpenBUGS.pdf) to
perform source attribution analysis on MLST data from *Campylobacter
jejuni*. This is so far just a mostly working sketch, if you try it, we
welcome feedback via issues here on github. 

## Installation

### Linux

Of course first you need to install R. In order to use this package
you also need to install OpenBUGS on your system. You can read about
how to do that [here](http://www.openbugs.net/w/Downloads) or you can
just run the script below in a terminal:

```sh
sudo apt-get install libc6-dev-i386
wget -O OpenBUGS.tar.gz "http://www.openbugs.net/w/OpenBUGS_3_2_3?action=AttachFile&do=get&target=OpenBUGS-3.2.3.tar.gz"
mkdir OpenBUGS
tar zxvf OpenBUGS.tar.gz -C OpenBUGS --strip-components=1
cd OpenBUGS
./configure
make
sudo make install
cd ..
```

Then install the package from github. You can clone it and install it
or you could use the following commands in R.

```r
library(ghit)
install_github("trosendal/kilde")
```

### Windows

  * Install R
  * Install OpenBUGS. Follow instructions
    [here](http://www.openbugs.net/w/Downloads). Or you can run the
    following powershell script to install it with all defaults:

```sh
Invoke-WebRequest "http://www.openbugs.net/w/OpenBUGS_3_2_3?action=AttachFile&do=get&target=OpenBUGS323setup.exe" -OutFile "..\OpenBUGS323setup.exe";Start-Process -FilePath "..\OpenBUGS323setup.exe" -ArgumentList "/VERYSILENT" -NoNewWindow -Wait
```

Then install the package from github. You can clone it and install it
or you could use the following commands in R to install it via the
`ghit` library. Either way you need to have the appropriate build
tools installed with R to install packages from source on
Windows. Read more about that
[here](https://cran.r-project.org/bin/windows/Rtools/)

```r
library(ghit)
install_github("trosendal/kilde")
```

If building the package doesn't work you should be able to install the
built binary. Download it from
[here](https://ci.appveyor.com/api/projects/trosendal/kilde/artifacts/kilde_0.1.zip). And
install it in R with the following commands:

```r
setwd("<The directory that you saved the file in>")
install.packages("kilde_0.1.zip", repos = NULL)
```

## Example

```r

##load the library
library(kilde)

## Get the sample data
df <- clean_data(country = "Peru")

## Construct inits for the analysis
inits <- function() {
    list(g0 = rep(1, df$bugs_data$ns), 
         pmuta = structure(.Data = rep(0.5, df$bugs_data$ns * 7), 
                           .Dim = c(df$bugs_data$ns, 7)), 
         h0=structure(.Data=rep(1, df$bugs_data$ns * df$bugs_data$ns), 
                      .Dim=c(df$bugs_data$ns, df$bugs_data$ns)))
					  }

## Define the parameters to monitor
parameters <- c("qASP", "phi", "etaASP")

## Run the default model
res <- run_model(df = df, 
                 inits = inits, 
                 parameters = parameters)

## Plot the result
par(mfrow = c(df$bugs_data$ns, 1),
    mar = c(5,5,1,1))
plotrows <- ceiling(df$bugs_data$ns ^ 0.5)
plotcols <- ceiling(df$bugs_data$ns / plotrows)
par(mfrow=c(plotrows, plotcols))
for(i in 1:df$bugs_data$ns){
    plot(density(res$sims.list$phi[,i]),
         main="",
         xlab=df$group_names[i],
         xlim=c(0,1))
    points(mean(res$sims.list$phi[,i]),
           0,
           col="red",
           pch=16)
}
```

## Wishlist

* Transition to JAGS
* Multiple well described models
* Remove as much BUGS code as possible ie. Do most of the data
  handling in R
