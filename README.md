[![Build Status](https://travis-ci.org/trosendal/kilde.svg)](https://travis-ci.org/trosendal/kilde)

# A package to perform Bayesian source attribution of Campylobacter in R

An R package that connects to
[OpenBUGS](http://www.openbugs.net) via
[R2OpenBUGS](https://cran.r-project.org/web/packages/R2OpenBUGS/vignettes/R2OpenBUGS.pdf) to
perform source attribution analysis on MLST data from *Campylobacter
jejuni*


## Installation

## Linux

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

Then install the package from github

```r
library(devtools)
install_github("trosendal/kilde")
```

## Example

```
library(kilde)
df <- sample_data()
```
