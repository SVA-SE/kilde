## Copyright (c) 2016, Jukka Ranta, Thomas Rosendal
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
## 1. Redistributions of source code must retain the above copyright
## notice, this list of conditions and the following disclaimer.
##
## 2. Redistributions in binary form must reproduce the above copyright
## notice, this list of conditions and the following disclaimer in the
## documentation and/or other materials provided with the distribution.
##
## 3. Neither the name of the copyright holder nor the names of its
## contributors may be used to endorse or promote products derived from
## this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
##' Read in the sample data and clean it
##' 
##' @title sample_data 
##' @return A data.frame
##' @author Thomas Rosendal
##' @import utils
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
