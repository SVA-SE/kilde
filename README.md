|Linux Build|Windows Build|Test Coverage|
|----|----|----|
|[![Build Status](https://travis-ci.org/trosendal/kilde.svg)](https://travis-ci.org/trosendal/kilde)|[![Build status](https://ci.appveyor.com/api/projects/status/dlxa7wxqs2lrc2eh?svg=true)](https://ci.appveyor.com/project/trosendal/kilde)|[![Coverage Status](https://coveralls.io/repos/github/trosendal/kilde/badge.svg?maxAge=600&branch=master)](https://coveralls.io/github/trosendal/kilde?branch=master)|
# A package to perform Bayesian source attribution of Campylobacter in R

An R package that connects to
[OpenBUGS](http://www.openbugs.net) via
[R2OpenBUGS](https://cran.r-project.org/web/packages/R2OpenBUGS/vignettes/R2OpenBUGS.pdf) to
perform source attribution analysis multiple character data. The
examples are focused on its application to MLST data from *Campylobacter
jejuni*. This is a working implementation of 2 models: One using just
the ST-type as the definition of genotype and the second using the 7
alleles, each as character definitions in the model. These two models
have both been implemented with the MCMC running in R and the MCMC
running in OpenBUGS.

The package contains a fake sample dataset that is used for testing
and generating examples. 

## Installation

The package is dependent on the software `OpenBUGS` and the R package
`R2OpenBUGS`, so you will need to install those before the package
works. The following sections 

### Windows

  * Install R
  * Install OpenBUGS. Follow instructions
    [here](http://www.openbugs.net/w/Downloads). Or you can run the
    following PowerShell script to install it with all defaults:

```sh
Invoke-WebRequest "http://www.openbugs.net/w/OpenBUGS_3_2_3?action=AttachFile&do=get&target=OpenBUGS323setup.exe" -OutFile "..\OpenBUGS323setup.exe";Start-Process -FilePath "..\OpenBUGS323setup.exe" -ArgumentList "/VERYSILENT" -NoNewWindow -Wait
```

Then install the package from github. When you install an R package,
you need to make sure you have 'write permission' to the library
location that you are installing to. If you do not have permission you
will get an error. The best way to avoid this problem is to use a
library location in your 'home' directory like:
`C:/Users/thomas.rosendal/Documents/`. You can read more about
managing libraries
[here](https://cran.r-project.org/doc/manuals/R-admin.html#Managing-libraries). To
actually install the package:

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
### Linux

Of course, first you need to install R. In order to use this package
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

### OSX

See the instructions [here](http://www.openbugs.net/w/Downloads) to
install OpenBUGS. After you have OpenBUGS you can clone and build the
package from source as explained in the Linux section above; however
this is untested.


## Usage

### Have a look at the Vignette

```r
library(kilde)
vignette("attribution", package = "kilde")
```

### An example of running the model with ST level character data

```r
## Load the library
library(kilde)
## Get some data
df <- sample_data()
## This dataset contains 4 different countries, we'll pick Canada:
df <- df[df$country == "Canada",]
## We need to format the data to the form accepted by the model function.
ob <- dataformatting_ST(DATA = df, UM = 2)
## Initialize the mcmc data objects. Note: 100 iterations is too few.
result <- initialize_mcmc_ST(ob$ns, ob$STu, MCMC = 100, ob$Nisolates)
## Run the MCMC for the ST model in R
mcmc_ob <- runmcmc_ST(result = result,
           ob = ob,
           h = 0,
           FULL = 0)
## plot the history of the MCMC
plot_history(mcmc_ob, 50)
## plot the modelfit
plot_modelfit(mcmc_ob, 50)
## A summary of the model fit
summary_kilde(mcmc_ob, 50)
## The sample attribution
plot_sample_attribution(mcmc_ob, 50)
## a plot of the population attribution
plot_population_attribution(mcmc_ob, 50)
##
##
## We can run the same ST model in BUGS
##
## Then the same model in BUGS
## Initialize the MCMC objects
initial_result <- initialize_bugs_ST(ob)
## Run the model, this time with 1000 iterations, because bugs will
## throw an error if we try to run too few.
result_bugs <- run_bugs(result = initial_result,
                        ob = ob,
                        MCMC = 1000,
                        n.burnin = 100,
                        FULL = 0,
                        model = "SA_ST_model.jag",
                        n.chains = 1)
plot_history(result_bugs, 100)
plot_modelfit(result_bugs, 100)
summary_kilde(result_bugs, 100)
plot_sample_attribution(result_bugs, 100)
plot_population_attribution(result_bugs, 100)
```
### The model estimated by using the Allele level results

```r
library(kilde)
## Read in a format data
df <- sample_data()
df <- df[df$country == "Canada",]
ob <- dataformatting(DATA = df,
                     UM = 2)
######################################
## Initialize and then run mcmc in R
######################################
result <- initialize_mcmc(ns = ob$inits$ns,
                          nat = ob$inits$nat,
                          MCMC = 100,
                          Nisolates = ob$inits$Nisolates)
mcmc_ob <- runmcmc(result, ob, MCMC = 100, h = 0, FULL = 0)
##  Plot the results of this model.
plot_history(mcmc_ob, 50)
plot_modelfit(mcmc_ob, 50)
summary_kilde(mcmc_ob, 50)
plot_sample_attribution(mcmc_ob, 50)
plot_population_attribution(mcmc_ob, 50)
##
################################################
## Initialize and then run the model in OpenBugs
################################################
##
## Below, BUGS model cannot handle a large number of MCMC iterations
## for all parameters.  Therefore, it is advisable to try with smaller
## number of iterations to start, perhaps 1000.
##
initial_result <- initialize_bugs(ob)
result_bugs <- run_bugs(result = initial_result,
                        ob = ob,
                        MCMC = 1000,
                        n.burnin = 100,
                        FULL = 0)
plot_history(result_bugs, 100)
plot_modelfit(result_bugs, 100)
summary_kilde(result_bugs, 100)
plot_sample_attribution(result_bugs, 100)
plot_population_attribution(result_bugs, 100)

```

## Wishlist

* Transition to JAGS
* Multiple well described models, including existing models like the
  island model and HALD model.
* Refactor code to reduce complexity of datacleaning/formatting
* Break large functions into smaller components
* Take more advantage of S3 classes to construct the internal wiring
  of the package
* When code has been substantially improved add a model for N-alleles

