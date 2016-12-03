rm(list=ls())

source('~/Dropbox/ChesBayTrends/code/wavk.test.new.R')
#source('~/Documents/Funtimes.Git.Rstudio/funtimes/funtimes/R/wavk.test.R')
V <- rnorm(100)

wavk.test.new.R(V~t, factor.length = "adaptive.selection", B=1000, out = TRUE)
