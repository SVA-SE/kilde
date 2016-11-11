|Linux Build|Windows Build|Test Coverage|
|----|----|----|
|[![Build Status](https://travis-ci.org/trosendal/kilde.svg)](https://travis-ci.org/trosendal/kilde)|[![Build status](https://ci.appveyor.com/api/projects/status/dlxa7wxqs2lrc2eh?svg=true)](https://ci.appveyor.com/project/trosendal/kilde)|[![Coverage Status](https://coveralls.io/repos/github/trosendal/kilde/badge.svg?maxAge=600&branch=master)](https://coveralls.io/github/trosendal/kilde?branch=master)|
# A package to perform Bayesian source attribution of Campylobacter in R

An R package that connects to
[OpenBUGS](http://www.openbugs.net) via
[R2OpenBUGS](https://cran.r-project.org/web/packages/R2OpenBUGS/vignettes/R2OpenBUGS.pdf) to
perform source attribution analysis on MLST data from *Campylobacter
jejuni*. This is so far just a mostly working sketch, if you try it, we
welcome feedback via issues here on github. 

## Installation

The package is dependent on the software `OpenBUGS` and the R package
`R2OpenBUGS`, so you will need to install those before the package
works. 

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

Then install the package from github. You can do this in a few
different ways depending on what you are comfortable with:

#### Use the `ghit` package to install from github source

Just run the following commands in R (You may not need to run the
first two lines: `install.packages()` if you already have those
packages installed)

```r
install.packages("R2OpenBUGS")
install.packages("ghit")
library(ghit)
install_github("trosendal/kilde")
```

#### Clone it build it, install it

1. Clone this repository with `git`
2. Navigate the the repository root and run `make install`
3. Install the R2OpenBUGS package in R


### Windows

  * Install R
  * Install OpenBUGS. Follow instructions
    [here](http://www.openbugs.net/w/Downloads). Or you can run the
    following PowerShell script to install it with all defaults:

```sh
Invoke-WebRequest "http://www.openbugs.net/w/OpenBUGS_3_2_3?action=AttachFile&do=get&target=OpenBUGS323setup.exe" -OutFile "..\OpenBUGS323setup.exe";Start-Process -FilePath "..\OpenBUGS323setup.exe" -ArgumentList "/VERYSILENT" -NoNewWindow -Wait
```

Then install the package from github. You can do this in a few
different ways depending on what you are comfortable with, choose one
of A, B or C:

#### A: Download prebuilt package:

1. Download the package from
[here](https://ci.appveyor.com/api/projects/trosendal/kilde/artifacts/kilde_0.1.zip)
2. Start R
3. Navigate the the directory that you saved the file to, perhaps your
   Downloads folder like this:

    ```r
    setwd("C:/Users/thomas.rosendal/Downloads")
    ```	
4. Install the package like this:

    ```r
    install.packages("kilde_0.1.zip", repos = NULL)
    ```

5. Install `R2OpenBUGS`

    ```r
    install.packages("R2OpenBUGS")
    ```

#### B: Use the `ghit` package to install from github source

To install from source you need the appropriate build tools for
Windows, you can read more about that
[here](https://cran.r-project.org/bin/windows/Rtools/). Once you have
the build tools installed you can run the following to install the package:

```r
install.packages("R2OpenBUGS")
library(ghit)
install_github("trosendal/kilde")
```

#### C: Clone it build it, install it

1. Clone this repository with `git`
2. Navigate the the repository root and run `make install`
3. Install the R2OpenBUGS package in R

## Usage

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
* Multiple well described models (This is currently being done...)
* Add tests (Work in progress)
* Remove as much BUGS code as possible ie. Do most of the data
  handling in R
