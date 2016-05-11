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
##' A method to run a model in bugs
##' @title run_model 
##' @param model The name of the model
##' This is the file name in the inst/models directory of the package
##' @param df the object from clean_data 
##' @param inits An inits function
##' @param parameters Parameters to monitor
##' @param n.chains The number of chains
##' @param n.burnin The burning length
##' @param n.iter The number of iterations
##' @return A model object
##' @author Thomas Rosendal
##' @import R2OpenBUGS
##' @export
run_model <- function(model = "SA_allelemodel1.jag",
                      df,
                      inits,
                      parameters,
                      n.chains = 1,
                      n.burnin = 99,
                      n.iter = 999) {
    bugs(df$bugs_data,
         inits,
         parameters,
         system.file(package = "kilde", paste0("models/", model)),
         n.chains = n.chains,
         n.burnin = n.burnin,
         n.iter = n.iter)
}
